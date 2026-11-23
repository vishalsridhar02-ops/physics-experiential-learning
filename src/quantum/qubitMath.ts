export type GateId = 'X' | 'Y' | 'Z' | 'H' | 'S' | 'I';
export type CardinalStateId = '0' | '1' | '+' | '-' | '+i' | '-i';

export interface Complex {
  re: number;
  im: number;
}

export interface QubitState {
  alpha: Complex;
  beta: Complex;
}

export interface CardinalDefinition {
  id: CardinalStateId;
  ketLabel: string;
  compactVec: string;
  description: string;
  state: QubitState;
}

export interface GateDefinition {
  id: GateId;
  name: string;
  role: string;
  scalarLabel: string | null;
  displayMatrix: [[string, string], [string, string]];
  matrix: [[Complex, Complex], [Complex, Complex]];
}

export interface MatrixStep {
  gateId: GateId;
  gate: GateDefinition;
  input: QubitState;
  output: QubitState;
  row0Terms: [string, string];
  row1Terms: [string, string];
  row0Simplified: [string, string];
  row1Simplified: [string, string];
  interpretation: string;
}

export interface MeasurementResult {
  shots: number;
  count0: number;
  count1: number;
  observedP0: number;
  observedP1: number;
  theoreticalP0: number;
  theoreticalP1: number;
  deltaP0: number;
  deltaP1: number;
  measuredKet: string;
}

export interface CircuitEvaluation {
  initialId: CardinalStateId;
  initialState: QubitState;
  steps: MatrixStep[];
  finalState: QubitState;
  previousState: QubitState | null;
}

const INV_SQRT2 = Math.SQRT1_2;
const HALF_INV_SQRT2 = Math.SQRT1_2 / 2;
const EPS = 1e-10;

const SNAP_VALUES = [
  0,
  1,
  -1,
  0.5,
  -0.5,
  INV_SQRT2,
  -INV_SQRT2,
  HALF_INV_SQRT2,
  -HALF_INV_SQRT2,
];

export function snapReal(v: number): number {
  for (let i = 0; i < SNAP_VALUES.length; i += 1) {
    if (Math.abs(v - SNAP_VALUES[i]) < EPS) {
      return SNAP_VALUES[i] === 0 ? 0 : SNAP_VALUES[i];
    }
  }
  return v;
}

export function complex(re: number, im = 0): Complex {
  return { re: snapReal(re), im: snapReal(im) };
}

export function cAdd(a: Complex, b: Complex): Complex {
  return complex(a.re + b.re, a.im + b.im);
}

export function cMul(a: Complex, b: Complex): Complex {
  return complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

export function cMagSq(z: Complex): number {
  const m2 = z.re * z.re + z.im * z.im;
  if (Math.abs(m2) < EPS) return 0;
  if (Math.abs(m2 - 0.5) < EPS) return 0.5;
  if (Math.abs(m2 - 1) < EPS) return 1;
  return m2;
}

export function cMag(z: Complex): number {
  return Math.sqrt(cMagSq(z));
}

export function cIsZero(z: Complex): boolean {
  return Math.abs(z.re) < EPS && Math.abs(z.im) < EPS;
}

export function cEquals(a: Complex, b: Complex): boolean {
  return Math.abs(a.re - b.re) < EPS && Math.abs(a.im - b.im) < EPS;
}

/**
 * Returns complex phase angle in [0, 360) degrees, or null if amplitude is 0.
 */
export function cPhaseDeg(z: Complex): number | null {
  if (cIsZero(z)) return null;
  let deg = (Math.atan2(z.im, z.re) * 180) / Math.PI;
  if (deg < -EPS) deg += 360;
  const rounded = Math.round(deg) % 360;
  if (Math.abs(deg - rounded) < 1e-6 || Math.abs(deg - 360) < 1e-6) {
    return rounded === 360 ? 0 : rounded;
  }
  return Number(deg.toFixed(1));
}

export function formatPhaseCardinal(deg: number | null): string {
  if (deg === null) return '—';
  if (deg === 0) return '0° (+1)';
  if (deg === 90) return '90° (+i)';
  if (deg === 180) return '180° (-1)';
  if (deg === 270) return '270° (-i)';
  return `${deg}°`;
}

export function getRelativePhaseDeg(state: QubitState): number | null {
  const pA = cPhaseDeg(state.alpha);
  const pB = cPhaseDeg(state.beta);
  if (pA === null || pB === null) return null;
  let diff = (pB - pA) % 360;
  if (diff < 0) diff += 360;
  if (Math.abs(diff - Math.round(diff)) < 1e-6) {
    return Math.round(diff) % 360;
  }
  return Number(diff.toFixed(1));
}

export interface RelativePhaseStatus {
  deltaPhi: number | null;
  label: string;
  isSuperposition: boolean;
}

export function getRelativePhaseStatus(state: QubitState): RelativePhaseStatus {
  if (cIsZero(state.alpha) || cIsZero(state.beta)) {
    return {
      deltaPhi: null,
      label: 'Basis state (no superposition)',
      isSuperposition: false,
    };
  }
  const rel = getRelativePhaseDeg(state);
  if (rel === 0) {
    return { deltaPhi: 0, label: 'In-phase superposition (Δφ = 0°)', isSuperposition: true };
  }
  if (rel === 90) {
    return { deltaPhi: 90, label: 'Positive imaginary phase (Δφ = 90° / +i)', isSuperposition: true };
  }
  if (rel === 180) {
    return { deltaPhi: 180, label: 'Phase-flipped superposition (Δφ = 180°)', isSuperposition: true };
  }
  if (rel === 270) {
    return { deltaPhi: 270, label: 'Negative imaginary phase (Δφ = 270° / -i)', isSuperposition: true };
  }
  return {
    deltaPhi: rel,
    label: `Superposition (Δφ = ${rel}°)`,
    isSuperposition: true,
  };
}

export function isPurelyRealState(state: QubitState): boolean {
  return Math.abs(state.alpha.im) < EPS && Math.abs(state.beta.im) < EPS;
}

const C0 = complex(0, 0);
const C1 = complex(1, 0);
const CM1 = complex(-1, 0);
const CI = complex(0, 1);
const CMI = complex(0, -1);
const CH = complex(INV_SQRT2, 0);
const CMH = complex(-INV_SQRT2, 0);
const CHI = complex(0, INV_SQRT2);
const CMHI = complex(0, -INV_SQRT2);

export const CARDINAL_ORDER: CardinalStateId[] = ['0', '1', '+', '-', '+i', '-i'];

export const CARDINAL_STATES: Record<CardinalStateId, CardinalDefinition> = {
  '0': {
    id: '0',
    ketLabel: '|0⟩',
    compactVec: '[1, 0]ᵀ',
    description: 'Computational basis |0⟩ (+Z)',
    state: { alpha: C1, beta: C0 },
  },
  '1': {
    id: '1',
    ketLabel: '|1⟩',
    compactVec: '[0, 1]ᵀ',
    description: 'Computational basis |1⟩ (-Z)',
    state: { alpha: C0, beta: C1 },
  },
  '+': {
    id: '+',
    ketLabel: '|+⟩',
    compactVec: '[1, 1]ᵀ/√2',
    description: 'Equal superposition |+⟩ (+X)',
    state: { alpha: CH, beta: CH },
  },
  '-': {
    id: '-',
    ketLabel: '|-⟩',
    compactVec: '[1, -1]ᵀ/√2',
    description: 'Phase-flipped superposition |-⟩ (-X)',
    state: { alpha: CH, beta: CMH },
  },
  '+i': {
    id: '+i',
    ketLabel: '|+i⟩',
    compactVec: '[1, i]ᵀ/√2',
    description: 'Complex superposition |+i⟩ (+Y)',
    state: { alpha: CH, beta: CHI },
  },
  '-i': {
    id: '-i',
    ketLabel: '|-i⟩',
    compactVec: '[1, -i]ᵀ/√2',
    description: 'Complex superposition |-i⟩ (-Y)',
    state: { alpha: CH, beta: CMHI },
  },
};

export const GATES: Record<GateId, GateDefinition> = {
  X: {
    id: 'X',
    name: 'Pauli-X',
    role: 'Bit-Flip',
    scalarLabel: null,
    displayMatrix: [
      ['0', '1'],
      ['1', '0'],
    ],
    matrix: [
      [C0, C1],
      [C1, C0],
    ],
  },
  Y: {
    id: 'Y',
    name: 'Pauli-Y',
    role: 'Bit + Phase',
    scalarLabel: null,
    displayMatrix: [
      ['0', '-i'],
      ['i', '0'],
    ],
    matrix: [
      [C0, CMI],
      [CI, C0],
    ],
  },
  Z: {
    id: 'Z',
    name: 'Pauli-Z',
    role: 'Phase-Flip',
    scalarLabel: null,
    displayMatrix: [
      ['1', '0'],
      ['0', '-1'],
    ],
    matrix: [
      [C1, C0],
      [C0, CM1],
    ],
  },
  H: {
    id: 'H',
    name: 'Hadamard',
    role: 'Superposition',
    scalarLabel: '1/√2',
    displayMatrix: [
      ['1', '1'],
      ['1', '-1'],
    ],
    matrix: [
      [CH, CH],
      [CH, CMH],
    ],
  },
  S: {
    id: 'S',
    name: 'Phase',
    role: 'π/2 Phase',
    scalarLabel: null,
    displayMatrix: [
      ['1', '0'],
      ['0', 'i'],
    ],
    matrix: [
      [C1, C0],
      [C0, CI],
    ],
  },
  I: {
    id: 'I',
    name: 'Identity',
    role: 'No-Op',
    scalarLabel: null,
    displayMatrix: [
      ['1', '0'],
      ['0', '1'],
    ],
    matrix: [
      [C1, C0],
      [C0, C1],
    ],
  },
};

export function createCardinalState(id: CardinalStateId): QubitState {
  const def = CARDINAL_STATES[id];
  return {
    alpha: { ...def.state.alpha },
    beta: { ...def.state.beta },
  };
}

function formatRealExact(v: number): string {
  const s = snapReal(v);
  if (s === 0) return '0';
  if (s === 1) return '1';
  if (s === -1) return '-1';
  if (s === 0.5) return '1/2';
  if (s === -0.5) return '-1/2';
  if (s === INV_SQRT2) return '1/√2';
  if (s === -INV_SQRT2) return '-1/√2';
  if (s === HALF_INV_SQRT2) return '1/(2√2)';
  if (s === -HALF_INV_SQRT2) return '-1/(2√2)';
  return s.toFixed(4);
}

function formatImagExact(v: number): string {
  const s = snapReal(v);
  if (s === 0) return '0';
  if (s === 1) return 'i';
  if (s === -1) return '-i';
  if (s === 0.5) return 'i/2';
  if (s === -0.5) return '-i/2';
  if (s === INV_SQRT2) return 'i/√2';
  if (s === -INV_SQRT2) return '-i/√2';
  if (s === HALF_INV_SQRT2) return 'i/(2√2)';
  if (s === -HALF_INV_SQRT2) return '-i/(2√2)';
  return `${s.toFixed(4)}i`;
}

export function formatComplexExact(z: Complex): string {
  const re = snapReal(z.re);
  const im = snapReal(z.im);
  if (re === 0 && im === 0) return '0';
  if (im === 0) return formatRealExact(re);
  if (re === 0) return formatImagExact(im);

  const absRe = Math.abs(re);
  const absIm = Math.abs(im);
  if (Math.abs(absRe - absIm) < EPS) {
    const rSign = re < 0 ? '-1' : '1';
    const iSign = im < 0 ? '-i' : '+i';
    if (absRe === 0.5) return `(${rSign}${iSign})/2`;
    if (absRe === INV_SQRT2) return `(${rSign}${iSign})/√2`;
    if (absRe === HALF_INV_SQRT2) return `(${rSign}${iSign})/(2√2)`;
    if (absRe === 1) return `${rSign}${iSign}`;
  }

  const rPart = formatRealExact(re);
  const iPart = formatImagExact(Math.abs(im));
  return `${rPart} ${im > 0 ? '+' : '-'} ${iPart}`;
}

export function formatComplexDecimal(z: Complex): string {
  const re = snapReal(z.re);
  const im = snapReal(z.im);
  if (re === 0 && im === 0) return '0';
  if (im === 0) {
    if (re === 1) return '1.0000';
    if (re === -1) return '-1.0000';
    return re.toFixed(4);
  }
  if (re === 0) {
    if (im === 1) return 'i';
    if (im === -1) return '-i';
    return `${im.toFixed(4)}i`;
  }
  const sign = im > 0 ? ' + ' : ' - ';
  return `${re.toFixed(4)}${sign}${Math.abs(im).toFixed(4)}i`;
}

export function getStateProbabilities(state: QubitState): {
  p0: number;
  p1: number;
  sum: number;
  p0Fraction: string;
  p1Fraction: string;
  p0Percent: string;
  p1Percent: string;
} {
  const p0 = cMagSq(state.alpha);
  const p1 = cMagSq(state.beta);
  const sum = p0 + p1;

  const toFraction = (p: number): string => {
    if (p === 0) return '0';
    if (p === 1) return '1';
    if (p === 0.5) return '1/2';
    return p.toFixed(4);
  };

  return {
    p0,
    p1,
    sum,
    p0Fraction: toFraction(p0),
    p1Fraction: toFraction(p1),
    p0Percent: `${(p0 * 100).toFixed(1)}%`,
    p1Percent: `${(p1 * 100).toFixed(1)}%`,
  };
}

function formatKetCoeff(z: Complex, isSecondTerm: boolean): { sign: string; body: string } {
  const exact = formatComplexExact(z);
  const re = snapReal(z.re);
  const im = snapReal(z.im);

  if (!isSecondTerm) {
    if (exact === '1') return { sign: '', body: '' };
    if (exact === '-1') return { sign: '-', body: '' };
    if (exact === 'i') return { sign: '', body: 'i' };
    if (exact === '-i') return { sign: '-', body: 'i' };
    if (re !== 0 && im !== 0) return { sign: '', body: `(${exact})` };
    if (exact.startsWith('-')) return { sign: '-', body: `(${exact.slice(1)})` };
    return { sign: '', body: `(${exact})` };
  }

  // Second term with binary + or - operator
  if (exact === '1') return { sign: ' + ', body: '' };
  if (exact === '-1') return { sign: ' - ', body: '' };
  if (exact === 'i') return { sign: ' + ', body: 'i' };
  if (exact === '-i') return { sign: ' - ', body: 'i' };
  if (re === 0 || im === 0) {
    if (exact.startsWith('-')) return { sign: ' - ', body: `(${exact.slice(1)})` };
    return { sign: ' + ', body: `(${exact})` };
  }
  return { sign: ' + ', body: `(${exact})` };
}

export function formatKetNotation(state: QubitState): {
  ket: string;
  alias: string | null;
  kind: string;
} {
  const { alpha, beta } = state;
  const aZero = cIsZero(alpha);
  const bZero = cIsZero(beta);

  let ket = '';
  if (aZero && !bZero) {
    const c = formatKetCoeff(beta, false);
    ket = `${c.sign}${c.body}|1⟩`;
  } else if (!aZero && bZero) {
    const c = formatKetCoeff(alpha, false);
    ket = `${c.sign}${c.body}|0⟩`;
  } else {
    const c0 = formatKetCoeff(alpha, false);
    const c1 = formatKetCoeff(beta, true);
    ket = `${c0.sign}${c0.body}|0⟩${c1.sign}${c1.body}|1⟩`;
  }

  // Identify exact cardinal alias or equivalence up to global phase
  if (aZero && !bZero) {
    const pB = cPhaseDeg(beta);
    if (pB === 0) return { ket, alias: null, kind: 'Computational basis |1⟩' };
    return { ket, alias: `≡ |1⟩ (phase ${formatPhaseCardinal(pB)})`, kind: 'Basis |1⟩ with global phase' };
  }
  if (!aZero && bZero) {
    const pA = cPhaseDeg(alpha);
    if (pA === 0) return { ket, alias: null, kind: 'Computational basis |0⟩' };
    return { ket, alias: `≡ |0⟩ (phase ${formatPhaseCardinal(pA)})`, kind: 'Basis |0⟩ with global phase' };
  }

  const p0 = cMagSq(alpha);
  const p1 = cMagSq(beta);
  const rel = getRelativePhaseDeg(state);
  const pA = cPhaseDeg(alpha);

  if (p0 === 0.5 && p1 === 0.5 && rel !== null) {
    const cardinalMap: Record<number, { name: string; factoredKet: string; label: string }> = {
      0: { name: '|+⟩', factoredKet: '(|0⟩ + |1⟩)/√2', label: 'In-phase superposition (Δφ = 0°)' },
      90: { name: '|+i⟩', factoredKet: '(|0⟩ + i|1⟩)/√2', label: 'Positive imaginary phase (Δφ = 90° / +i)' },
      180: { name: '|-⟩', factoredKet: '(|0⟩ - |1⟩)/√2', label: 'Phase-flipped superposition (Δφ = 180°)' },
      270: { name: '|-i⟩', factoredKet: '(|0⟩ - i|1⟩)/√2', label: 'Negative imaginary phase (Δφ = 270° / -i)' },
    };
    const match = cardinalMap[rel];
    if (match) {
      if (pA === 0) {
        return { ket: match.factoredKet, alias: `= ${match.name}`, kind: match.label };
      }
      return {
        ket,
        alias: `≡ ${match.name} (global phase ${pA}°)`,
        kind: match.label,
      };
    }
  }

  return { ket, alias: null, kind: 'Superposition state' };
}

export function getStateAngleDegrees(state: QubitState): number {
  if (isPurelyRealState(state)) {
    const deg = (Math.atan2(state.beta.re, state.alpha.re) * 180) / Math.PI;
    const rounded = Math.round(deg);
    return Math.abs(deg - rounded) < 1e-9 ? (rounded === 0 ? 0 : rounded) : Number(deg.toFixed(1));
  }
  const deg = (Math.atan2(cMag(state.beta), cMag(state.alpha)) * 180) / Math.PI;
  const rounded = Math.round(deg);
  return Math.abs(deg - rounded) < 1e-9 ? (rounded === 0 ? 0 : rounded) : Number(deg.toFixed(1));
}

function describeStep(gateId: GateId, input: QubitState, output: QubitState): string {
  const inInfo = formatKetNotation(input);
  const outInfo = formatKetNotation(output);
  const inLbl = inInfo.alias ? inInfo.alias.replace(/^[=≡]\s*/, '').split(' ')[0] : inInfo.ket;
  const outLbl = outInfo.alias ? `${outInfo.ket} (${outInfo.alias.replace(/^[=≡]\s*/, '')})` : outInfo.ket;

  if (gateId === 'I') {
    return `Identity (I) leaves ${inLbl} invariant with unchanged magnitudes and phases.`;
  }
  if (gateId === 'X') {
    return `Pauli-X exchanges the |0⟩ and |1⟩ amplitudes, transforming ${inLbl} into ${outLbl}.`;
  }
  if (gateId === 'Y') {
    return `Pauli-Y swaps |0⟩ and |1⟩ while applying ±i phase shifts, mapping ${inLbl} to ${outLbl}.`;
  }
  if (gateId === 'Z') {
    return `Pauli-Z preserves |0⟩ and applies a 180° (-1) phase shift to |1⟩, transforming ${inLbl} into ${outLbl}.`;
  }
  if (gateId === 'S') {
    return `Phase gate (S) preserves |0⟩ and rotates the |1⟩ phase by +90° (×i), transforming ${inLbl} into ${outLbl}.`;
  }
  return `Hadamard (H) interferes α and β with (α+β)/√2 and (α-β)/√2, mapping ${inLbl} to ${outLbl}.`;
}

export function applyGate(state: QubitState, gateId: GateId): MatrixStep {
  const gate = GATES[gateId];
  const [[m00, m01], [m10, m11]] = gate.matrix;

  const p00 = cMul(m00, state.alpha);
  const p01 = cMul(m01, state.beta);
  const p10 = cMul(m10, state.alpha);
  const p11 = cMul(m11, state.beta);

  const output: QubitState = {
    alpha: cAdd(p00, p01),
    beta: cAdd(p10, p11),
  };

  const aStr = formatComplexExact(state.alpha);
  const bStr = formatComplexExact(state.beta);

  const row0Terms: [string, string] = [
    `(${formatComplexExact(m00)})(${aStr})`,
    `(${formatComplexExact(m01)})(${bStr})`,
  ];
  const row1Terms: [string, string] = [
    `(${formatComplexExact(m10)})(${aStr})`,
    `(${formatComplexExact(m11)})(${bStr})`,
  ];

  const row0Simplified: [string, string] = [formatComplexExact(p00), formatComplexExact(p01)];
  const row1Simplified: [string, string] = [formatComplexExact(p10), formatComplexExact(p11)];

  return {
    gateId,
    gate,
    input: state,
    output,
    row0Terms,
    row1Terms,
    row0Simplified,
    row1Simplified,
    interpretation: describeStep(gateId, state, output),
  };
}

export function evaluateCircuit(
  initialId: CardinalStateId,
  gateSequence: GateId[]
): CircuitEvaluation {
  const initialState = createCardinalState(initialId);
  const steps: MatrixStep[] = [];
  let cursor = initialState;

  for (let i = 0; i < gateSequence.length; i += 1) {
    const step = applyGate(cursor, gateSequence[i]);
    steps.push(step);
    cursor = step.output;
  }

  return {
    initialId,
    initialState,
    steps,
    finalState: cursor,
    previousState: steps.length > 0 ? steps[steps.length - 1].input : null,
  };
}

export function simulateMeasurements(
  state: QubitState,
  shots = 100,
  rng: () => number = Math.random
): MeasurementResult {
  const { p0, p1 } = getStateProbabilities(state);
  let count0 = 0;
  let count1 = 0;

  for (let i = 0; i < shots; i += 1) {
    if (rng() < p0) {
      count0 += 1;
    } else {
      count1 += 1;
    }
  }

  const observedP0 = count0 / shots;
  const observedP1 = count1 / shots;

  return {
    shots,
    count0,
    count1,
    observedP0,
    observedP1,
    theoreticalP0: p0,
    theoreticalP1: p1,
    deltaP0: observedP0 - p0,
    deltaP1: observedP1 - p1,
    measuredKet: formatKetNotation(state).ket,
  };
}
