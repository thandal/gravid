#!/usr/bin/env node

"use strict";

const crypto = require("node:crypto");
const fs = require("node:fs");
const Module = require("node:module");
const os = require("node:os");
const path = require("node:path");
const { pathToFileURL } = require("node:url");

const EXTENSION_PREFIX = "yzane.markdown-pdf-";

function usage() {
  console.log(`Usage: node scripts/build_article_pdf.cjs [SOURCE.md] [OUTPUT.pdf]

Build a Markdown manuscript with the VS Code Markdown PDF renderer.

Defaults:
  SOURCE.md   article.md
  OUTPUT.pdf  SOURCE with its extension replaced by .pdf

Optional environment variables:
  MARKDOWN_PDF_EXTENSION  Path to an installed yzane.markdown-pdf extension
  MARKDOWN_PDF_CHROME     Path to Chrome, Chromium, or Edge`);
}

function fail(message) {
  throw new Error(message);
}

function compareVersions(left, right) {
  const parse = (value) => value.split(/[.-]/).map((part) => {
    const number = Number.parseInt(part, 10);
    return Number.isNaN(number) ? part : number;
  });
  const a = parse(left);
  const b = parse(right);
  const length = Math.max(a.length, b.length);
  for (let index = 0; index < length; index += 1) {
    const av = a[index] ?? 0;
    const bv = b[index] ?? 0;
    if (av === bv) continue;
    if (typeof av === "number" && typeof bv === "number") return av - bv;
    return String(av).localeCompare(String(bv));
  }
  return 0;
}

function extensionVersion(extensionPath) {
  try {
    const manifest = JSON.parse(
      fs.readFileSync(path.join(extensionPath, "package.json"), "utf8"),
    );
    return manifest.version || "0";
  } catch {
    return "0";
  }
}

function findMarkdownPdfExtension() {
  if (process.env.MARKDOWN_PDF_EXTENSION) {
    const explicit = path.resolve(process.env.MARKDOWN_PDF_EXTENSION);
    if (!fs.existsSync(path.join(explicit, "dist", "extension.js"))) {
      fail(
        `MARKDOWN_PDF_EXTENSION does not contain dist/extension.js: ${explicit}`,
      );
    }
    return explicit;
  }

  const extensionRoots = [
    path.join(os.homedir(), ".vscode-server", "extensions"),
    path.join(os.homedir(), ".vscode", "extensions"),
    path.join(os.homedir(), ".cursor-server", "extensions"),
  ];
  const candidates = [];
  for (const root of extensionRoots) {
    if (!fs.existsSync(root)) continue;
    for (const entry of fs.readdirSync(root, { withFileTypes: true })) {
      if (!entry.isDirectory() || !entry.name.startsWith(EXTENSION_PREFIX)) {
        continue;
      }
      const candidate = path.join(root, entry.name);
      if (fs.existsSync(path.join(candidate, "dist", "extension.js"))) {
        candidates.push(candidate);
      }
    }
  }

  if (candidates.length === 0) {
    fail(
      "Could not find the VS Code Markdown PDF extension (yzane.markdown-pdf). " +
        "Install it or set MARKDOWN_PDF_EXTENSION to its installation directory.",
    );
  }
  candidates.sort((left, right) =>
    compareVersions(extensionVersion(right), extensionVersion(left)),
  );
  return candidates[0];
}

function setNested(target, keys, value) {
  let cursor = target;
  for (let index = 0; index < keys.length - 1; index += 1) {
    const key = keys[index];
    if (!cursor[key] || typeof cursor[key] !== "object") cursor[key] = {};
    cursor = cursor[key];
  }
  cursor[keys[keys.length - 1]] = value;
}

function buildConfiguration(manifest, stagingDirectory) {
  const sections = {
    "markdown-pdf": {},
    markdown: { styles: [] },
    http: { proxy: "" },
  };
  const properties = manifest.contributes?.configuration?.properties || {};
  const prefix = "markdown-pdf.";
  for (const [fullKey, specification] of Object.entries(properties)) {
    if (!fullKey.startsWith(prefix)) continue;
    setNested(
      sections["markdown-pdf"],
      fullKey.slice(prefix.length).split("."),
      specification.default,
    );
  }
  sections["markdown-pdf"].outputDirectory = stagingDirectory;
  if (process.env.MARKDOWN_PDF_CHROME) {
    const chrome = path.resolve(process.env.MARKDOWN_PDF_CHROME);
    if (!fs.existsSync(chrome)) {
      fail(`MARKDOWN_PDF_CHROME does not exist: ${chrome}`);
    }
    sections["markdown-pdf"].executablePath = chrome;
  }
  return sections;
}

class Uri {
  constructor(scheme, fsPath, raw) {
    this.scheme = scheme;
    this.fsPath = fsPath;
    this.path = fsPath;
    this.raw = raw;
  }

  static file(filePath) {
    const absolute = path.resolve(filePath);
    return new Uri("file", absolute, pathToFileURL(absolute).href);
  }

  static parse(value) {
    if (/^[A-Za-z][A-Za-z0-9+.-]*:/.test(value)) {
      const parsed = new URL(value);
      const fsPath = parsed.protocol === "file:"
        ? decodeURIComponent(parsed.pathname)
        : parsed.pathname;
      return new Uri(parsed.protocol.slice(0, -1), fsPath, value);
    }
    return new Uri("", value, value);
  }

  toString() {
    return this.raw;
  }
}

function createVscodeShim(source, workspaceRoot, sections, registered, errors) {
  const disposable = () => ({ dispose() {} });
  const configuration = (section) => {
    const values = sections[section] || {};
    return Object.assign({}, values, {
      get(key, fallback) {
        return values[key] === undefined ? fallback : values[key];
      },
    });
  };
  const log = (label, message) => {
    process.stderr.write(`[markdown-pdf${label}] ${message}\n`);
  };
  const outputChannel = {
    info(message) { log("", message); },
    warn(message) { log(" warning", message); },
    error(message) { log(" error", message); },
    show() {},
    dispose() {},
  };

  return {
    version: "build-script",
    Uri,
    ProgressLocation: { Notification: 15 },
    env: {
      language: "en",
      openExternal: async () => true,
    },
    commands: {
      registerCommand(name, callback) {
        registered.set(name, callback);
        return disposable();
      },
      executeCommand: async () => undefined,
    },
    workspace: {
      getConfiguration: configuration,
      getWorkspaceFolder() {
        return { uri: Uri.file(workspaceRoot) };
      },
      onDidSaveTextDocument() {
        return disposable();
      },
    },
    window: {
      activeTextEditor: {
        document: {
          languageId: "markdown",
          uri: Uri.file(source),
          fileName: source,
          isUntitled: false,
          getText() {
            return fs.readFileSync(source, "utf8");
          },
        },
      },
      createOutputChannel() {
        return outputChannel;
      },
      setStatusBarMessage() {
        return disposable();
      },
      showInformationMessage: async (message) => log("", message),
      showWarningMessage: async (message) => log(" warning", message),
      showErrorMessage: async (message) => {
        errors.push(String(message));
        log(" error", message);
      },
      withProgress: async (_options, task) => task(),
    },
  };
}

function validatePdf(pdfPath) {
  const stat = fs.statSync(pdfPath);
  if (stat.size < 10_000) {
    fail(`Generated PDF is unexpectedly small (${stat.size} bytes): ${pdfPath}`);
  }
  const descriptor = fs.openSync(pdfPath, "r");
  try {
    const signature = Buffer.alloc(5);
    fs.readSync(descriptor, signature, 0, signature.length, 0);
    if (signature.toString("ascii") !== "%PDF-") {
      fail(`Generated file does not have a PDF signature: ${pdfPath}`);
    }
  } finally {
    fs.closeSync(descriptor);
  }
  return stat;
}

async function build(source, output) {
  if (!fs.existsSync(source)) fail(`Source file does not exist: ${source}`);
  if (path.extname(source).toLowerCase() !== ".md") {
    fail(`Source file must have a .md extension: ${source}`);
  }
  if (path.extname(output).toLowerCase() !== ".pdf") {
    fail(`Output file must have a .pdf extension: ${output}`);
  }

  const extensionRoot = findMarkdownPdfExtension();
  const manifest = JSON.parse(
    fs.readFileSync(path.join(extensionRoot, "package.json"), "utf8"),
  );
  const stagingDirectory = fs.mkdtempSync(
    path.join(os.tmpdir(), "gravid-article-pdf-"),
  );
  const outputDirectory = path.dirname(output);
  const temporaryOutput = path.join(
    outputDirectory,
    `.${path.basename(output)}.${process.pid}.tmp`,
  );

  try {
    const sections = buildConfiguration(manifest, stagingDirectory);
    const registered = new Map();
    const errors = [];
    const vscode = createVscodeShim(
      source,
      path.dirname(source),
      sections,
      registered,
      errors,
    );

    const originalLoad = Module._load;
    Module._load = function patchedLoad(request, parent, isMain) {
      if (request === "vscode") return vscode;
      return originalLoad.call(this, request, parent, isMain);
    };

    let extension;
    try {
      extension = require(path.join(extensionRoot, "dist", "extension.js"));
    } finally {
      Module._load = originalLoad;
    }

    const context = {
      subscriptions: [],
      extension: { packageJSON: manifest },
      globalStorageUri: Uri.file(
        path.join(os.tmpdir(), "gravid-markdown-pdf-cache"),
      ),
    };
    extension.activate(context);
    const render = registered.get("extension.markdown-pdf.pdf");
    if (!render) fail("Markdown PDF export command was not registered");
    await render();
    if (errors.length > 0) fail(errors.join("\n"));

    const stagedPdf = path.join(
      stagingDirectory,
      `${path.basename(source, path.extname(source))}.pdf`,
    );
    const stat = validatePdf(stagedPdf);
    fs.mkdirSync(outputDirectory, { recursive: true });
    fs.copyFileSync(stagedPdf, temporaryOutput);
    fs.renameSync(temporaryOutput, output);

    const digest = crypto
      .createHash("sha256")
      .update(fs.readFileSync(output))
      .digest("hex");
    console.log(
      `Built ${path.relative(process.cwd(), output) || output} ` +
        `(${stat.size} bytes, SHA-256 ${digest})`,
    );
    console.log(
      `Renderer: Markdown PDF ${manifest.version} from ${extensionRoot}`,
    );
  } finally {
    if (fs.existsSync(temporaryOutput)) fs.rmSync(temporaryOutput, { force: true });
    fs.rmSync(stagingDirectory, { recursive: true, force: true });
  }
}

if (process.argv.includes("--help") || process.argv.includes("-h")) {
  usage();
  process.exit(0);
}

const source = path.resolve(process.argv[2] || "article.md");
const defaultOutput = source.replace(/\.[^.]+$/, ".pdf");
const output = path.resolve(process.argv[3] || defaultOutput);

build(source, output).catch((error) => {
  process.stderr.write(`PDF build failed: ${error.message}\n`);
  process.exitCode = 1;
});
