# Original repository and public release provenance

## Repository locations

- Access-controlled original repository: <https://github.com/lifrary-01/tpd-navigator>
- Curated public repository: <https://github.com/Mr-BBAWWI/tpd-navigator-submission-2026>
- Public demonstration URL: <https://mr-bbawwi.github.io/tpd-navigator-submission-2026/>

The original repository remains private, with its existing history preserved under owner-controlled access; separately authorized and reviewed changes may be published to an authorized existing branch in that repository. The public repository is the destination for a new history containing a curated source snapshot; this exporter itself does not change the original source, history, Git data, or remotes, copy `.git`, invoke Git, contact a remote service, or deploy itself. Public demonstration hosting and Git publication of the curated code are separate release actions.

## Curated boundary

The public snapshot includes allowlisted application, package, contract, test and script source; selected public design/reference structures; vendor license/provenance material and the approved Vina executable; and the curated saved static demonstration. It excludes private expert evidence, original papers and office documents, credentials, local databases, logs, caches, reports, request/response captures, private builders, model weights, raw ensembles, archives and historical repository metadata.

Full private evidence packages and historical material require explicit access from the repository owner. Their absence must not be interpreted as synthetic replacement evidence.

## Authorized export-only privacy transformations

The exporter never edits the original files. In the exported copy only, it:

1. replaces an authorized reviewer display name with `Reviewer A` in exactly three reviewed source files while preserving identifiers, authorization logic, computation and approval behavior; and
2. removes the `expert_documents` metadata key from the public copy of the design-source manifest because that key describes excluded private documents.

For each transformed file, `source-export-manifest.json` records the original source SHA-256, exported SHA-256 and a category-level transformation reason. It does not list excluded private filenames. All other allowlisted payloads are copied byte-for-byte.

## Verification limits

The export manifest closes over regular payload files but excludes itself from its own hash closure. Verification detects missing, additional, modified or linked files and repeats privacy/security scanning. This provenance statement documents packaging boundaries; it does not claim complete scientific approval, reproduce private frozen evidence, or verify public hosting availability.
