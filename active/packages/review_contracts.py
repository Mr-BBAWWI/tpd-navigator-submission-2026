"""Runtime bridge to the existing, independently versioned I3 review contracts."""
import importlib.util
from packages.contracts import ROOT, encoded

spec = importlib.util.spec_from_file_location("tpd_i3_contracts", ROOT / "contracts/review/v0.1.0/validation.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ReviewReader, validate, URI = module.ReviewReader, module.validate, module.URI


def save(port, kind, context, **fields):
    value = {"payload_type": kind, "payload_version": "0.1.0", **context, **fields}
    validate(kind, value)
    provenance = "test_fixture" if context["data_mode"] == "test_fixture" else "computed"
    return port.put_raw(encoded(value), "application/json", provenance, URI + "#/$defs/" + kind)
