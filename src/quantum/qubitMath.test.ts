import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import {
  CARDINAL_ORDER,
  applyGate,
  cMag,
  cPhaseDeg,
  createCardinalState,
  evaluateCircuit,
  formatComplexDecimal,
  formatComplexExact,
  formatKetNotation,
  getRelativePhaseDeg,
  getRelativePhaseStatus,
  getStateAngleDegrees,
  getStateProbabilities,
  simulateMeasurements,
} from './qubitMath';

test('All 6 cardinal initial states are normalized with exact expected phases', () => {
  for (const id of CARDINAL_ORDER) {
    const s = createCardinalState(id);
    const probs = getStateProbabilities(s);
    assert.equal(probs.sum.toFixed(4), '1.0000', `Cardinal state ${id} must be normalized`);
  }

  assert.equal(getRelativePhaseDeg(createCardinalState('+')), 0);
  assert.equal(getRelativePhaseDeg(createCardinalState('-')), 180);
  assert.equal(getRelativePhaseDeg(createCardinalState('+i')), 90);
  assert.equal(getRelativePhaseDeg(createCardinalState('-i')), 270);

  assert.equal(formatComplexExact(createCardinalState('+i').beta), 'i/√2');
  assert.equal(formatComplexDecimal(createCardinalState('+i').beta), '0.7071i');
  assert.equal(formatComplexExact(createCardinalState('-i').beta), '-i/√2');
  assert.equal(formatComplexDecimal(createCardinalState('-i').beta), '-0.7071i');

  // Clean factored symbolic ket notation for the 4 canonical superpositions
  assert.equal(formatKetNotation(createCardinalState('+')).ket, '(|0⟩ + |1⟩)/√2');
  assert.equal(formatKetNotation(createCardinalState('-')).ket, '(|0⟩ - |1⟩)/√2');
  assert.equal(formatKetNotation(createCardinalState('+i')).ket, '(|0⟩ + i|1⟩)/√2');
  assert.equal(formatKetNotation(createCardinalState('-i')).ket, '(|0⟩ - i|1⟩)/√2');

  // Relative Phase Status exact descriptive labels + magnitudes
  assert.equal(getRelativePhaseStatus(createCardinalState('0')).label, 'Basis state (no superposition)');
  assert.equal(getRelativePhaseStatus(createCardinalState('1')).label, 'Basis state (no superposition)');
  assert.equal(getRelativePhaseStatus(createCardinalState('+')).label, 'In-phase superposition (Δφ = 0°)');
  assert.equal(getRelativePhaseStatus(createCardinalState('+i')).label, 'Positive imaginary phase (Δφ = 90° / +i)');
  assert.equal(getRelativePhaseStatus(createCardinalState('-')).label, 'Phase-flipped superposition (Δφ = 180°)');
  assert.equal(getRelativePhaseStatus(createCardinalState('-i')).label, 'Negative imaginary phase (Δφ = 270° / -i)');
  assert.equal(cMag(createCardinalState('+i').alpha).toFixed(4), '0.7071');
  assert.equal(cMag(createCardinalState('+i').beta).toFixed(4), '0.7071');
});

test('Pauli-X, Pauli-Y, Pauli-Z, Hadamard, Phase S, and Identity gates satisfy exact complex identities', () => {
  // X|0> = |1>, X|1> = |0>
  assert.equal(formatKetNotation(evaluateCircuit('0', ['X']).finalState).ket, '|1⟩');
  assert.equal(formatKetNotation(evaluateCircuit('1', ['X']).finalState).ket, '|0⟩');

  // Y|0> = i|1>, Y|1> = -i|0>
  const y0 = evaluateCircuit('0', ['Y']).finalState;
  assert.equal(formatKetNotation(y0).ket, 'i|1⟩');
  assert.equal(formatComplexExact(y0.beta), 'i');
  assert.equal(formatComplexDecimal(y0.beta), 'i');
  assert.equal(cPhaseDeg(y0.beta), 90);
  assert.equal(getStateProbabilities(y0).p1, 1);

  const y1 = evaluateCircuit('1', ['Y']).finalState;
  assert.equal(formatKetNotation(y1).ket, '-i|0⟩');
  assert.equal(formatComplexExact(y1.alpha), '-i');
  assert.equal(formatComplexDecimal(y1.alpha), '-i');
  assert.equal(cPhaseDeg(y1.alpha), 270);
  assert.equal(getStateProbabilities(y1).p0, 1);

  // Y|+i> = |+i> (eigenstate with +1 eigenvalue)
  const yPlusI = evaluateCircuit('+i', ['Y']).finalState;
  assert.equal(formatKetNotation(yPlusI).alias, '= |+i⟩');
  assert.equal(getRelativePhaseDeg(yPlusI), 90);

  // Z|+> = |-> and Z|+i> = |-i>
  const zPlus = evaluateCircuit('+', ['Z']).finalState;
  assert.equal(formatKetNotation(zPlus).alias, '= |-⟩');
  assert.equal(cPhaseDeg(zPlus.beta), 180);

  const zPlusI = evaluateCircuit('+i', ['Z']).finalState;
  assert.equal(formatKetNotation(zPlusI).alias, '= |-i⟩');
  assert.equal(cPhaseDeg(zPlusI.beta), 270);

  // S|+> = |+i>, S(S|+>) = |-> (S^2 = Z), S|+i> = |->, S|-> = |-i>
  const sPlus = evaluateCircuit('+', ['S']).finalState;
  assert.equal(formatKetNotation(sPlus).alias, '= |+i⟩');
  assert.equal(cPhaseDeg(sPlus.beta), 90);

  const ssPlus = evaluateCircuit('+', ['S', 'S']).finalState;
  assert.equal(formatKetNotation(ssPlus).alias, '= |-⟩');
  assert.equal(cPhaseDeg(ssPlus.beta), 180);

  const sMinus = evaluateCircuit('-', ['S']).finalState;
  assert.equal(formatKetNotation(sMinus).alias, '= |-i⟩');
  assert.equal(cPhaseDeg(sMinus.beta), 270);

  // H(H|0>) = |0> and H|+i> -> ((1+i)/2)|0> + ((1-i)/2)|1> (relative phase 270° = -90°, normalized)
  const hh0 = evaluateCircuit('0', ['H', 'H']).finalState;
  assert.equal(formatKetNotation(hh0).ket, '|0⟩');
  assert.equal(getStateAngleDegrees(hh0), 0);

  const hPlusI = evaluateCircuit('+i', ['H']).finalState;
  assert.equal(formatComplexExact(hPlusI.alpha), '(1+i)/2');
  assert.equal(formatComplexExact(hPlusI.beta), '(1-i)/2');
  assert.equal(getRelativePhaseDeg(hPlusI), 270);
  const probsHPlusI = getStateProbabilities(hPlusI);
  assert.equal(probsHPlusI.p0, 0.5);
  assert.equal(probsHPlusI.p1, 0.5);
  assert.equal(probsHPlusI.sum.toFixed(4), '1.0000');

  // Identity leaves complex states untouched
  const iPlusI = applyGate(createCardinalState('+i'), 'I').output;
  assert.equal(formatKetNotation(iPlusI).ket, formatKetNotation(createCardinalState('+i')).ket);
});

test('Complex decimal formatting never outputs messy 0.0000 + 1.0000i strings', () => {
  const y0 = evaluateCircuit('0', ['Y']).finalState;
  assert.equal(formatComplexDecimal(y0.alpha), '0');
  assert.equal(formatComplexDecimal(y0.beta), 'i');

  const h1 = evaluateCircuit('1', ['H']).finalState;
  assert.equal(formatComplexDecimal(h1.alpha), '0.7071');
  assert.equal(formatComplexDecimal(h1.beta), '-0.7071');

  const hs0 = evaluateCircuit('0', ['H', 'S']).finalState;
  assert.equal(formatComplexDecimal(hs0.alpha), '0.7071');
  assert.equal(formatComplexDecimal(hs0.beta), '0.7071i');
});

test('100-shot measurement simulator uses complex squared magnitudes accurately', () => {
  const y0 = evaluateCircuit('0', ['Y']).finalState;
  const mY0 = simulateMeasurements(y0, 100);
  assert.equal(mY0.count0, 0);
  assert.equal(mY0.count1, 100);

  let tick = 0;
  const detRng = () => (tick++ % 2 === 0 ? 0.25 : 0.75);
  const plusI = createCardinalState('+i');
  const mPlusI = simulateMeasurements(plusI, 100, detRng);
  assert.equal(mPlusI.count0, 50);
  assert.equal(mPlusI.count1, 50);
  assert.equal(mPlusI.deltaP0, 0);
});

test('Self-contained index.html inline vanilla JS complex kernel matches all Clifford identities', () => {
  const htmlPath = path.resolve(process.cwd(), 'index.html');
  const html = fs.readFileSync(htmlPath, 'utf-8');
  const scriptMatch = html.match(/<script>([\s\S]*?)<\/script>/);
  assert.ok(scriptMatch, 'index.html should contain an inline <script> block');
  assert.equal(/<script[^>]+src=/.test(html), false, 'index.html must be self-contained');

  const sandbox: { window: Record<string, unknown>; Math: Math } = {
    window: {},
    Math,
  };
  vm.createContext(sandbox);
  vm.runInContext(scriptMatch[1], sandbox);

  const kernel = sandbox.window.__QUBIT_KERNEL__ as {
    CARDINAL_ORDER: string[];
    evaluateCircuit: (
      init: string,
      seq: string[]
    ) => { finalState: { alpha: { re: number; im: number }; beta: { re: number; im: number } }; steps: unknown[] };
    formatKetNotation: (s: unknown) => { ket: string; alias: string | null };
    formatComplexDecimal: (z: unknown) => string;
    cMag: (z: unknown) => number;
    cPhaseDeg: (z: unknown) => number | null;
    getRelativePhaseDeg: (s: unknown) => number | null;
    getRelativePhaseStatus: (s: unknown) => { deltaPhi: number | null; label: string; isSuperposition: boolean };
    getStateProbabilities: (s: unknown) => { p0: number; p1: number; sum: number };
    simulateMeasurements: (s: unknown, n: number) => { count0: number; count1: number };
  };

  assert.ok(kernel, 'inline script should expose __QUBIT_KERNEL__');
  assert.deepEqual(Array.from(kernel.CARDINAL_ORDER), ['0', '1', '+', '-', '+i', '-i']);

  // H then S on |0> -> |+i> with factored ket, arg(beta) = 90 deg, and exact Relative Phase Status
  const hs = kernel.evaluateCircuit('0', ['H', 'S']);
  assert.equal(kernel.formatKetNotation(hs.finalState).ket, '(|0⟩ + i|1⟩)/√2');
  assert.equal(kernel.formatKetNotation(hs.finalState).alias, '= |+i⟩');
  assert.equal(kernel.formatComplexDecimal(hs.finalState.beta), '0.7071i');
  assert.equal(kernel.cMag(hs.finalState.beta).toFixed(4), '0.7071');
  assert.equal(kernel.cPhaseDeg(hs.finalState.alpha), 0);
  assert.equal(kernel.cPhaseDeg(hs.finalState.beta), 90);
  assert.equal(kernel.getRelativePhaseStatus(hs.finalState).label, 'Positive imaginary phase (Δφ = 90° / +i)');
  assert.equal(kernel.getStateProbabilities(hs.finalState).sum.toFixed(4), '1.0000');

  // Y|0> = i|1> -> basis state (no superposition)
  const y0 = kernel.evaluateCircuit('0', ['Y']);
  assert.equal(kernel.formatKetNotation(y0.finalState).ket, 'i|1⟩');
  assert.equal(kernel.cPhaseDeg(y0.finalState.beta), 90);
  assert.equal(kernel.getRelativePhaseStatus(y0.finalState).label, 'Basis state (no superposition)');

  // Z|+i> = |-i> -> Negative imaginary phase (Δφ = 270° / -i)
  const zpi = kernel.evaluateCircuit('+i', ['Z']);
  assert.equal(kernel.formatKetNotation(zpi.finalState).ket, '(|0⟩ - i|1⟩)/√2');
  assert.equal(kernel.formatKetNotation(zpi.finalState).alias, '= |-i⟩');
  assert.equal(kernel.getRelativePhaseDeg(zpi.finalState), 270);
  assert.equal(kernel.getRelativePhaseStatus(zpi.finalState).label, 'Negative imaginary phase (Δφ = 270° / -i)');

  // H|1> = |-> -> Phase-flipped superposition (Δφ = 180°)
  const h1 = kernel.evaluateCircuit('1', ['H']);
  assert.equal(kernel.formatKetNotation(h1.finalState).ket, '(|0⟩ - |1⟩)/√2');
  assert.equal(kernel.getRelativePhaseStatus(h1.finalState).label, 'Phase-flipped superposition (Δφ = 180°)');

  // H|+i> normalized
  const hpi = kernel.evaluateCircuit('+i', ['H']);
  assert.equal(kernel.getStateProbabilities(hpi.finalState).sum.toFixed(4), '1.0000');
  assert.equal(kernel.getRelativePhaseDeg(hpi.finalState), 270);
  assert.equal(kernel.getRelativePhaseStatus(hpi.finalState).label, 'Negative imaginary phase (Δφ = 270° / -i)');

  // 100 shots sum to 100
  const shots = kernel.simulateMeasurements(hs.finalState, 100);
  assert.equal(shots.count0 + shots.count1, 100);
});
