export type GateId = 'X' | 'H' | 'I';
export type BasisState = '0' | '1';

/**
 * Exact algebraic representation for single-qubit real stabilizer amplitudes:
 * value = (r / 2) + (s / (2 * sqrt(2)))
 * Under {X, H, I} starting from |0> or |1>, all states remain in:
 * 0 -> (0,0), ±1 -> (±2,0), ±1/√2 -> (0,±2)
 */
export interface ExactAmp {
  r: number;
  s: number;
}

export interface QubitState {
  alpha: ExactAmp;
  beta: ExactAmp;
}

export interface GateDefinition {
  id: GateId;
  name: string;
  shortLabel: string;
  role: string;
  scalarLabel: string | null;
  displayMatrix: [[string, string], [string, string]];
  exactMatrix: [[ExactAmp, ExactAmp], [ExactAmp, ExactAmp]];
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

const INV_SQRT2 = Math.SQRT1_2;

export const ZERO: ExactAmp = { r: 0, s: 0 };
export const ONE: ExactAmp = { r: 2, s: 0 };
export const NEG_ONE: ExactAmp = { r: -2, s: 0 };
export const INV_SQRT2_EXACT: ExactAmp = { r: 0, s: 2 };
export const NEG_INV_SQRT2_EXACT: ExactAmp = { r: 0, s: -2 };

export const GATES: Record<GateId, GateDefinition> = {
  X: {
    id: 'X',
    name: 'Pauli-X',
    shortLabel: 'X',
    role: 'Bit-Flip',
    scalarLabel: null,
    displayMatrix: [
      ['0', '1'],
      ['1', '0'],
    ],
    exactMatrix: [
      [ZERO, ONE],
      [ONE, ZERO],
    ],
  },
  H: {
    id: 'H',
    name: 'Hadamard',
    shortLabel: 'H',
    role: 'Superposition',
    scalarLabel: '1/√2',
    displayMatrix: [
      ['1', '1'],
      ['1', '-1'],
    ],
    exactMatrix: [
      [INV_SQRT2_EXACT, INV_SQRT2_EXACT],
      [INV_SQRT2_EXACT, NEG_INV_SQRT2_EXACT],
    ],
  },
  I: {
    id: 'I',
    name: 'Identity',
    shortLabel: 'I',
    role: 'No-Op',
    scalarLabel: null,
    displayMatrix: [
      ['1', '0'],
      ['0', '1'],
    ],
    exactMatrix: [
      [ONE, ZERO],
      [ZERO, ONE],
    ],
  },
};

export function createBasisState(basis: BasisState): QubitState {
  return basis === '0'
    ? { alpha: { ...ONE }, beta: { ...ZERO } }
    : { alpha: { ...ZERO }, beta: { ...ONE } };
}

export function ampToFloat(a: ExactAmp): number {
  return a.r / 2 + (a.s / 2) * INV_SQRT2;
}

export function ampEquals(a: ExactAmp, b: ExactAmp): boolean {
  return a.r === b.r && a.s === b.s;
}

export function mulExact(a: ExactAmp, b: ExactAmp): ExactAmp {
  // (a.r/2 + a.s/(2√2)) * (b.r/2 + b.s/(2√2))
  // = (2*a.r*b.r + a.s*b.s) / 8 + (a.r*b.s + a.s*b.r) / (4√2)
  // For our closed set, numerator divides cleanly into halves.
  const rNum = 2 * a.r * b.r + a.s * b.s;
  const sNum = a.r * b.s + a.s * b.r;
  return {
    r: rNum / 4,
    s: sNum / 2,
  };
}

export function addExact(a: ExactAmp, b: ExactAmp): ExactAmp {
  return {
    r: a.r + b.r,
    s: a.s + b.s,
  };
}

export function formatExactAmp(a: ExactAmp): string {
  if (a.r === 0 && a.s === 0) return '0';
  if (a.s === 0) {
    if (a.r === 2) return '1';
    if (a.r === -2) return '-1';
    if (a.r === 1) return '1/2';
    if (a.r === -1) return '-1/2';
    if (a.r % 2 === 0) return String(a.r / 2);
    return `${a.r}/2`;
  }
  if (a.r === 0) {
    if (a.s === 2) return '1/√2';
    if (a.s === -2) return '-1/√2';
    if (a.s === 1) return '1/(2√2)';
    if (a.s === -1) return '-1/(2√2)';
    return `${a.s}/(2√2)`;
  }
  return ampToFloat(a).toFixed(4);
}

export function formatSignedDecimal(a: ExactAmp): string {
  const val = ampToFloat(a);
  if (Math.abs(val) < 1e-12) return '0.0000';
  const sign = val > 0 ? '+' : '-';
  return `${sign}${Math.abs(val).toFixed(4)}`;
}

export function formatUnsignedDecimal(a: ExactAmp): string {
  const val = ampToFloat(a);
  if (Math.abs(val) < 1e-12) return '0.0000';
  return val.toFixed(4);
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
  const a = ampToFloat(state.alpha);
  const b = ampToFloat(state.beta);
  // Exact squared magnitudes for our stabilizer states
  const p0 = Math.abs(a * a) < 1e-12 ? 0 : Math.abs(a * a - 0.5) < 1e-12 ? 0.5 : Math.abs(a * a - 1) < 1e-12 ? 1 : a * a;
  const p1 = Math.abs(b * b) < 1e-12 ? 0 : Math.abs(b * b - 0.5) < 1e-12 ? 0.5 : Math.abs(b * b - 1) < 1e-12 ? 1 : b * b;
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

export function formatKetNotation(state: QubitState): {
  ket: string;
  alias: string | null;
} {
  const { alpha, beta } = state;

  if (ampEquals(alpha, ONE) && ampEquals(beta, ZERO)) {
    return { ket: '|0⟩', alias: null };
  }
  if (ampEquals(alpha, NEG_ONE) && ampEquals(beta, ZERO)) {
    return { ket: '-|0⟩', alias: 'Global phase π' };
  }
  if (ampEquals(alpha, ZERO) && ampEquals(beta, ONE)) {
    return { ket: '|1⟩', alias: null };
  }
  if (ampEquals(alpha, ZERO) && ampEquals(beta, NEG_ONE)) {
    return { ket: '-|1⟩', alias: 'Global phase π' };
  }
  if (ampEquals(alpha, INV_SQRT2_EXACT) && ampEquals(beta, INV_SQRT2_EXACT)) {
    return { ket: '(1/√2)|0⟩ + (1/√2)|1⟩', alias: '|+⟩' };
  }
  if (ampEquals(alpha, INV_SQRT2_EXACT) && ampEquals(beta, NEG_INV_SQRT2_EXACT)) {
    return { ket: '(1/√2)|0⟩ - (1/√2)|1⟩', alias: '|-⟩' };
  }
  if (ampEquals(alpha, NEG_INV_SQRT2_EXACT) && ampEquals(beta, INV_SQRT2_EXACT)) {
    return { ket: '-(1/√2)|0⟩ + (1/√2)|1⟩', alias: '-|-⟩' };
  }
  if (ampEquals(alpha, NEG_INV_SQRT2_EXACT) && ampEquals(beta, NEG_INV_SQRT2_EXACT)) {
    return { ket: '-(1/√2)|0⟩ - (1/√2)|1⟩', alias: '-|+⟩' };
  }

  const aStr = formatExactAmp(alpha);
  const bFloat = ampToFloat(beta);
  const bSign = bFloat >= 0 ? ' + ' : ' - ';
  const bMag = formatExactAmp({ r: Math.abs(beta.r), s: Math.abs(beta.s) });
  return { ket: `(${aStr})|0⟩${bSign}(${bMag})|1⟩`, alias: null };
}

export function getStateAngleDegrees(state: QubitState): number {
  const a = ampToFloat(state.alpha);
  const b = ampToFloat(state.beta);
  const deg = (Math.atan2(b, a) * 180) / Math.PI;
  const rounded = Math.round(deg);
  return Math.abs(deg - rounded) < 1e-9 ? (rounded === -0 ? 0 : rounded) : Number(deg.toFixed(1));
}

function describeStep(gateId: GateId, input: QubitState, output: QubitState): string {
  const inKet = formatKetNotation(input).alias || formatKetNotation(input).ket;
  const outInfo = formatKetNotation(output);
  const outKet = outInfo.alias ? `${outInfo.ket} (${outInfo.alias})` : outInfo.ket;

  if (gateId === 'I') {
    return `Identity leaves ${inKet} unchanged with identical amplitudes and probabilities.`;
  }
  if (gateId === 'X') {
    return `Pauli-X swaps the |0⟩ and |1⟩ amplitudes (reflection across +45°), mapping ${inKet} to ${outKet}.`;
  }
  // Hadamard
  const inIsBasis = ampEquals(input.beta, ZERO) || ampEquals(input.alpha, ZERO);
  if (inIsBasis) {
    const angle = getStateAngleDegrees(output);
    const signStr = angle > 0 ? `+${angle}°` : `${angle}°`;
    return `Hadamard rotates basis state ${inKet} into equal superposition ${outKet} at ${signStr}.`;
  }
  return `Second Hadamard causes constructive/destructive interference, collapsing superposition ${inKet} back to basis state ${outKet}.`;
}

export function applyGate(state: QubitState, gateId: GateId): MatrixStep {
  const gate = GATES[gateId];
  const [[m00, m01], [m10, m11]] = gate.exactMatrix;

  const p00 = mulExact(m00, state.alpha);
  const p01 = mulExact(m01, state.beta);
  const p10 = mulExact(m10, state.alpha);
  const p11 = mulExact(m11, state.beta);

  const outAlpha = addExact(p00, p01);
  const outBeta = addExact(p10, p11);
  const output: QubitState = { alpha: outAlpha, beta: outBeta };

  const aStr = formatExactAmp(state.alpha);
  const bStr = formatExactAmp(state.beta);

  const row0Terms: [string, string] = [
    `(${formatExactAmp(m00)})(${aStr})`,
    `(${formatExactAmp(m01)})(${bStr})`,
  ];
  const row1Terms: [string, string] = [
    `(${formatExactAmp(m10)})(${aStr})`,
    `(${formatExactAmp(m11)})(${bStr})`,
  ];

  const row0Simplified: [string, string] = [formatExactAmp(p00), formatExactAmp(p01)];
  const row1Simplified: [string, string] = [formatExactAmp(p10), formatExactAmp(p11)];

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

export interface CircuitEvaluation {
  initialBasis: BasisState;
  initialState: QubitState;
  steps: MatrixStep[];
  finalState: QubitState;
  previousState: QubitState | null;
}

export function evaluateCircuit(
  initialBasis: BasisState,
  gateSequence: GateId[]
): CircuitEvaluation {
  const initialState = createBasisState(initialBasis);
  const steps: MatrixStep[] = [];
  let cursor = initialState;

  for (let i = 0; i < gateSequence.length; i += 1) {
    const step = applyGate(cursor, gateSequence[i]);
    steps.push(step);
    cursor = step.output;
  }

  return {
    initialBasis,
    initialState,
    steps,
    finalState: cursor,
    previousState: steps.length > 0 ? steps[steps.length - 1].input : null,
  };
}

