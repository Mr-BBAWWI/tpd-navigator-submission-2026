const OPS = new Set(['candidate_compare', 'failure_analysis', 'acceptance_plan']);
const E3S = new Set(['CRBN', 'VHL']);
const DANGEROUS = new Set(['__proto__', 'prototype', 'constructor', 'url', 'command']);

export class ValidationError extends Error {
  constructor(code, message) {
    super(message);
    this.code = code;
  }
}

function plain(value) {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return false;
  const p = Object.getPrototypeOf(value);
  return p === Object.prototype || p === null;
}

function safeKeys(value, allowed) {
  if (!plain(value)) throw new ValidationError('INVALID_REQUEST', '요청 객체 형식이 올바르지 않습니다.');
  for (const key of Object.keys(value)) {
    if (DANGEROUS.has(key) || !allowed.has(key)) throw new ValidationError('UNKNOWN_FIELD', '허용되지 않은 요청 필드가 있습니다.');
  }
}

function nonEmptyString(value) {
  return typeof value === 'string' && value.length > 0 && value.trim() === value;
}

function finiteNumber(value) {
  return typeof value === 'number' && Number.isFinite(value);
}

function finiteInt(value) {
  return Number.isSafeInteger(value) && value >= 0;
}

export function canonicalStringify(value) {
  const seen = new Set();
  function visit(v) {
    if (v === null || typeof v === 'boolean' || typeof v === 'string') return JSON.stringify(v);
    if (typeof v === 'number') {
      if (!Number.isFinite(v)) throw new TypeError('Non-finite number');
      return JSON.stringify(v);
    }
    if (Array.isArray(v)) return '[' + v.map(visit).join(',') + ']';
    if (typeof v !== 'object') throw new TypeError('Unsupported canonical value');
    if (seen.has(v)) throw new TypeError('Cyclic value');
    seen.add(v);
    const keys = Object.keys(v).sort();
    const out = '{' + keys.map(k => JSON.stringify(k) + ':' + visit(v[k])).join(',') + '}';
    seen.delete(v);
    return out;
  }
  return visit(value);
}

const PARENT_COUNT_FIELDS = ['analog_records', 'assemblies', 'computed_but_pose_filter_rejected', 'docking_attempts', 'docking_completed', 'docking_skipped', 'failed_docking', 'selected_analogs'];

function holdoutOf(receipt) {
  const h = receipt.known_crbn_holdout ?? receipt.evidence.known_crbn_holdout;
  if (!plain(h)) throw new TypeError('Invalid holdout');
  for (const field of ['selected_pass_count', 'baseline_pass_count', 'raw_pass_count', 'raw_count']) {
    if (!finiteInt(h[field])) throw new TypeError('Invalid holdout count');
  }
  if (!Array.isArray(h.selected_per_seed_metrics) || h.selected_per_seed_metrics.length === 0) throw new TypeError('Invalid holdout rows');
  if (h.selected_pass_count > h.selected_per_seed_metrics.length || h.baseline_pass_count > h.selected_per_seed_metrics.length || h.raw_pass_count > h.raw_count) throw new TypeError('Invalid holdout count relationship');
  const visit = value => {
    if (value === null || typeof value === 'string' || typeof value === 'boolean') return;
    if (typeof value === 'number') {
      if (!Number.isFinite(value)) throw new TypeError('Invalid holdout value');
      return;
    }
    if (Array.isArray(value)) {
      for (const item of value) visit(item);
      return;
    }
    if (!plain(value)) throw new TypeError('Invalid holdout value');
    for (const [key, item] of Object.entries(value)) {
      if (/pass(?:es|ed)?(?:_all)?(?:_predicates)?$/i.test(key) && typeof item !== 'boolean') throw new TypeError('Invalid holdout predicate');
      visit(item);
    }
  };
  for (const row of h.selected_per_seed_metrics) {
    if (!plain(row)) throw new TypeError('Invalid holdout row');
    visit(row);
  }
  return h;
}

export function validateReceipt(receipt) {
  if (!plain(receipt) || !plain(receipt.public_scope) || !plain(receipt.evidence) || receipt.public_scope.scientific_approved !== false || receipt.public_scope.human_review_performed !== false || receipt.evidence.scientific_accepted !== false) throw new TypeError('Invalid frozen receipt');
  const e = receipt.evidence;
  if (!Array.isArray(e.parents) || !Array.isArray(e.candidates) || !Array.isArray(e.analogs)) throw new TypeError('Invalid evidence arrays');
  const snap = e.assessment_snapshot;
  if (!plain(snap) || snap.original_scientific_accepted !== false || !Array.isArray(snap.criteria14) || snap.criteria14.length !== 14) throw new TypeError('Invalid assessment snapshot');
  const criterionIDs = new Set();
  for (const criterion of snap.criteria14) {
    if (!plain(criterion) || !nonEmptyString(criterion.id) || criterionIDs.has(criterion.id) || !new Set(['pass', 'failed', 'pending']).has(criterion.status)) throw new TypeError('Invalid assessment criterion');
    criterionIDs.add(criterion.id);
  }
  const parentIDs = new Set();
  for (const p of e.parents) {
    if (!plain(p) || !nonEmptyString(p.parent_id) || parentIDs.has(p.parent_id) || !plain(p.counts_recomputed_from_records)) throw new TypeError('Invalid parent row');
    for (const field of PARENT_COUNT_FIELDS) if (!finiteInt(p.counts_recomputed_from_records[field])) throw new TypeError('Invalid parent count');
    parentIDs.add(p.parent_id);
  }
  const candidateKeys = new Set();
  for (const c of e.candidates) {
    if (!plain(c) || !nonEmptyString(c.candidate_id) || !nonEmptyString(c.candidate_key) || candidateKeys.has(c.candidate_key) || !nonEmptyString(c.parent_id) || !parentIDs.has(c.parent_id) || !E3S.has(c.e3_type) || !nonEmptyString(c.linker_id)) throw new TypeError('Invalid candidate row');
    candidateKeys.add(c.candidate_key);
  }
  holdoutOf(receipt);
  return true;
}

export function validateRequest(body, receipt) {
  validateReceipt(receipt);
  safeKeys(body, new Set(['operation', 'filters']));
  if (!OPS.has(body.operation)) throw new ValidationError('UNKNOWN_OPERATION', '지원하지 않는 작업입니다.');
  const filters = body.filters === undefined ? {} : body.filters;
  safeKeys(filters, new Set(['parent_id', 'e3_type', 'candidate_ids']));
  if (filters.parent_id !== undefined && !nonEmptyString(filters.parent_id)) throw new ValidationError('INVALID_FILTER', 'parent_id 형식이 올바르지 않습니다.');
  if (filters.e3_type !== undefined && !E3S.has(filters.e3_type)) throw new ValidationError('INVALID_FILTER', 'e3_type 형식이 올바르지 않습니다.');
  if (filters.candidate_ids !== undefined) {
    if (!Array.isArray(filters.candidate_ids) || filters.candidate_ids.length < 1 || filters.candidate_ids.length > 3) throw new ValidationError('INVALID_FILTER', 'candidate_ids는 1개 이상 3개 이하이어야 합니다.');
    const ids = new Set();
    for (const id of filters.candidate_ids) {
      if (!nonEmptyString(id)) throw new ValidationError('INVALID_FILTER', 'candidate_id 형식이 올바르지 않습니다.');
      if (ids.has(id)) throw new ValidationError('DUPLICATE_CANDIDATE_ID', 'candidate_id가 중복되었습니다.');
      ids.add(id);
    }
  }
  const candidates = receipt.evidence.candidates;
  const parentIDs = new Set(receipt.evidence.parents.map(p => p.parent_id));
  if (filters.parent_id !== undefined && !parentIDs.has(filters.parent_id)) throw new ValidationError('UNKNOWN_PARENT_ID', 'receipt에 없는 parent_id입니다.');
  if (filters.candidate_ids) {
    for (const id of filters.candidate_ids) {
      const all = candidates.filter(c => c.candidate_id === id);
      if (all.length === 0) throw new ValidationError('UNKNOWN_CANDIDATE_ID', 'receipt에 없는 candidate_id입니다.');
      const scoped = all.filter(c => (filters.parent_id === undefined || c.parent_id === filters.parent_id) && (filters.e3_type === undefined || c.e3_type === filters.e3_type));
      if (scoped.length > 1 || (scoped.length === 0 && all.length > 1 && filters.parent_id === undefined && filters.e3_type === undefined) || (all.length > 1 && filters.parent_id === undefined && filters.e3_type === undefined)) {
        throw new ValidationError('AMBIGUOUS_CANDIDATE_ID', 'candidate_id가 여러 receipt 레코드와 일치하므로 parent_id 또는 e3_type이 필요합니다.');
      }
    }
  }
  return Object.freeze({
    operation: body.operation,
    filters: Object.freeze({
      ...(filters.parent_id === undefined ? {} : { parent_id: filters.parent_id }),
      ...(filters.e3_type === undefined ? {} : { e3_type: filters.e3_type }),
      ...(filters.candidate_ids === undefined ? {} : { candidate_ids: Object.freeze([...filters.candidate_ids]) })
    })
  });
}

function candidateScope(receipt, filters) {
  let rows = receipt.evidence.candidates.slice();
  if (filters.parent_id !== undefined) rows = rows.filter(c => c.parent_id === filters.parent_id);
  if (filters.e3_type !== undefined) rows = rows.filter(c => c.e3_type === filters.e3_type);
  if (filters.candidate_ids !== undefined) {
    const ids = new Set(filters.candidate_ids);
    rows = rows.filter(c => ids.has(c.candidate_id));
  }
  return rows.sort((a, b) => a.candidate_key.localeCompare(b.candidate_key));
}

function numberOrNull(value) {
  return finiteNumber(value) ? value : null;
}

function publicCandidate(c) {
  const m = plain(c.selected_analog_docking_metrics) ? c.selected_analog_docking_metrics : {};
  return {
    evidence_id: 'E_CANDIDATES',
    candidate_key: c.candidate_key,
    candidate_id: c.candidate_id,
    parent_id: c.parent_id,
    e3_type: c.e3_type,
    linker_id: typeof c.linker_id === 'string' ? c.linker_id : null,
    orientation: typeof c.orientation === 'string' ? c.orientation : null,
    risk_flags: Array.isArray(c.risk_flags) ? c.risk_flags.filter(x => typeof x === 'string') : [],
    canonical_smiles: typeof c.canonical_smiles === 'string' ? c.canonical_smiles : null,
    selected: c.selected === true,
    strict_eligible: c.strict_eligible === true,
    pose_qualified: c.pose_qualified === true,
    pose_metrics: {
      best_core_rmsd_A: numberOrNull(m.best_core_rmsd_A),
      passing_pose_count: finiteInt(m.passing_pose_count) ? m.passing_pose_count : 0,
      pose_count: finiteInt(m.pose_count) ? m.pose_count : 0
    }
  };
}

function chooseCandidates(rows) {
  if (rows.length <= 3) return rows;
  const byParent = new Map();
  for (const row of rows) {
    if (!byParent.has(row.parent_id)) byParent.set(row.parent_id, { CRBN: [], VHL: [] });
    byParent.get(row.parent_id)[row.e3_type].push(row);
  }
  for (const parent of [...byParent.keys()].sort()) {
    const group = byParent.get(parent);
    if (group.CRBN.length && group.VHL.length) {
      const selected = [group.CRBN[0], group.VHL[0]];
      const used = new Set(selected.map(x => x.candidate_key));
      const third = rows.find(x => !used.has(x.candidate_key));
      if (third) selected.push(third);
      return selected.sort((a, b) => a.candidate_key.localeCompare(b.candidate_key));
    }
  }
  return rows.slice(0, 3);
}

function compare(receipt, filters) {
  const rows = candidateScope(receipt, filters);
  const chosen = chooseCandidates(rows);
  const rmsds = rows.map(c => c.selected_analog_docking_metrics?.best_core_rmsd_A).filter(finiteNumber);
  const byE3 = { CRBN: 0, VHL: 0 };
  const byParent = {};
  const linkers = new Set();
  for (const c of rows) {
    byE3[c.e3_type]++;
    byParent[c.parent_id] = (byParent[c.parent_id] || 0) + 1;
    if (typeof c.linker_id === 'string') linkers.add(c.linker_id);
  }
  return {
    evidence_ids: ['E_COUNTS', 'E_CANDIDATES'],
    count_scope: 'candidate records matching the supplied filters; filters are never broadened',
    no_match: rows.length === 0,
    no_match_explanation: rows.length === 0 ? 'Known filters are contradictory or select no saved candidate records.' : null,
    counts: {
      matched_candidates: rows.length,
      distinct_parents: Object.keys(byParent).length,
      parent_counts: byParent,
      e3_counts: byE3,
      distinct_linkers: linkers.size,
      selected: rows.filter(c => c.selected === true).length,
      strict_eligible: rows.filter(c => c.strict_eligible === true).length,
      pose_qualified: rows.filter(c => c.pose_qualified === true).length
    },
    metric_summary: {
      best_core_rmsd_A_min: rmsds.length ? Math.min(...rmsds) : null,
      best_core_rmsd_A_max: rmsds.length ? Math.max(...rmsds) : null,
      interpretation_limit: 'Saved per-candidate pose metric range; it is not a cross-E3 efficacy or ranking score.'
    },
    displayed_candidates: chosen.map(publicCandidate),
    selection_rule: 'candidate_key lexical order, preferring one CRBN and one VHL from the same parent when available; display selection is not a scientific ranking'
  };
}

function sumParentCounts(parents) {
  const fields = ['analog_records', 'assemblies', 'computed_but_pose_filter_rejected', 'docking_attempts', 'docking_completed', 'docking_skipped', 'failed_docking', 'selected_analogs'];
  const out = Object.fromEntries(fields.map(f => [f, 0]));
  for (const p of parents) {
    const counts = p.counts_recomputed_from_records;
    for (const f of fields) {
      const n = counts[f];
      if (!finiteInt(n)) throw new TypeError('Invalid parent count');
      out[f] += n;
    }
  }
  return out;
}

function failures(receipt, filters) {
  const candidates = candidateScope(receipt, filters);
  let parentIDs;
  if (filters.parent_id !== undefined && filters.e3_type === undefined && filters.candidate_ids === undefined) parentIDs = new Set([filters.parent_id]);
  else if (filters.e3_type !== undefined || filters.candidate_ids !== undefined) parentIDs = new Set(candidates.map(c => c.parent_id));
  else parentIDs = new Set(receipt.evidence.parents.map(p => p.parent_id));
  const parents = receipt.evidence.parents.filter(p => parentIDs.has(p.parent_id)).sort((a, b) => a.parent_id.localeCompare(b.parent_id));
  const samples = [];
  for (const a of receipt.evidence.analogs) if (parentIDs.has(a.parent_id) && typeof a.docking_failure === 'string' && a.docking_failure) samples.push({ record_id: a.analog_id, parent_id: a.parent_id, reason: a.docking_failure });
  for (const c of receipt.evidence.candidates) if (parentIDs.has(c.parent_id)) {
    if (typeof c.assembly_failure === 'string' && c.assembly_failure) samples.push({ record_id: c.candidate_id, parent_id: c.parent_id, reason: c.assembly_failure });
    if (typeof c.selected_analog_docking_failure === 'string' && c.selected_analog_docking_failure) samples.push({ record_id: c.candidate_id, parent_id: c.parent_id, reason: c.selected_analog_docking_failure });
  }
  samples.sort((a, b) => (a.parent_id + '/' + a.record_id + '/' + a.reason).localeCompare(b.parent_id + '/' + b.record_id + '/' + b.reason));
  const counts = sumParentCounts(parents);
  return {
    evidence_ids: ['E_COUNTS', 'E_FAILURES'],
    count_scope: 'parent-level saved docking counters; E3/candidate selection selects parent set only',
    selected_parent_ids: parents.map(p => p.parent_id),
    selected_candidate_count: candidates.length,
    counts,
    invariants: {
      completed_plus_failed_relationship_is_not_inferred: true,
      failed_docking: counts.failed_docking,
      computed_but_pose_filter_rejected: counts.computed_but_pose_filter_rejected
    },
    reason_availability: {
      nonnull_saved_reason_count: samples.length,
      missing_causes_are_not_inferred: true
    },
    failure_reason_samples: samples.slice(0, 5),
    limitation: 'Saved parent counters cannot attribute docking failures to an E3 branch, and null failure strings are not converted into invented causes.'
  };
}

function compact(value, max = 700) {
  const text = canonicalStringify(value);
  return text.length <= max ? value : { summary: text.slice(0, max) + '…', truncated: true };
}

const NEXT = {
  parent_funnel: 'more_candidate_breadth', expert_parent_selection: 'trusted_expert_parent_selection',
  actual_broad_families: 'more_candidate_breadth', qualified_panel: 'more_candidate_breadth',
  core_interaction_preservation: 'core_microstate_quantitative_validation', microstates_h_direction: 'core_microstate_quantitative_validation',
  known_crbn_calibration: 'independent_calibration', novel_ternary_repeats: 'novel_geometry_review',
  novel_ternary_geometry: 'novel_geometry_review', exact_synthesis_review: 'exact_graph_synthesis_review',
  formal_expert_decision: 'formal_expert_review'
};

function acceptance(receipt, filters) {
  const snap = receipt.evidence.assessment_snapshot;
  const rows = snap.criteria14.map(c => ({
    evidence_id: 'E_CRITERIA', id: c.id, title: c.title, status: c.status,
    reason: typeof c.reason === 'string' ? c.reason.slice(0, 240) : '', required: compact(c.required)
  }));
  const status_counts = { pass: 0, failed: 0, pending: 0 };
  for (const row of rows) {
    if (!(row.status in status_counts)) throw new TypeError('Invalid criterion status');
    status_counts[row.status]++;
  }
  const blockers = rows.filter(r => r.status !== 'pass').map(r => ({ id: r.id, status: r.status, title: r.title }));
  const next_steps = blockers.map(r => ({ criterion_id: r.id, action: NEXT[r.id] || 'criterion_specific_evidence_review', effect: 'collect or review evidence only; does not create a pass or approval' }));
  const h = holdoutOf(receipt);
  return {
    evidence_ids: ['E_CRITERIA', 'E_HOLDOUT', 'E_CANDIDATES'],
    immutable_snapshot: { assessment_id: snap.assessment_id, revision: snap.revision, original_scientific_accepted: snap.original_scientific_accepted === true },
    status_counts,
    criteria: rows,
    blockers,
    next_steps,
    selected_scope: { candidate_count: candidateScope(receipt, filters).length, filters },
    holdout_summary: {
      evidence_id: 'E_HOLDOUT',
      selected_pass_count: h.selected_pass_count,
      selected_total: h.selected_per_seed_metrics.length,
      baseline_pass_count: h.baseline_pass_count,
      raw_pass_count: h.raw_pass_count,
      raw_count: h.raw_count,
      scope_limit: 'Limited holdout on the same 6BOY development structure; it does not override the formal 14-criterion snapshot.'
    },
    scientific_approval: false,
    planning_limit: 'Planning steps are newly computed suggestions and are not executed science, status writes, or gate promotion.'
  };
}

export function computeEvidence(receipt, normalized) {
  validateReceipt(receipt);
  if (!plain(normalized) || !OPS.has(normalized.operation) || !plain(normalized.filters)) throw new TypeError('Invalid normalized request');
  if (normalized.operation === 'candidate_compare') return compare(receipt, normalized.filters);
  if (normalized.operation === 'failure_analysis') return failures(receipt, normalized.filters);
  return acceptance(receipt, normalized.filters);
}
