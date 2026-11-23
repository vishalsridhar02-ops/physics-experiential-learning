import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import {
  applyGate,
  createBasisState,
  evaluateCircuit,
  formatExactAmp,
  formatKetNotation,
  formatSignedDecimal,
  getStateAngleDegrees,
  getStateProbabilities,
  simulateMeasurements,
} from './qubitMath';

test('Pauli-X flips |0> to |1> and |1> to |0>', () => {
  const s0 = createBasisState('0');
  const stepX0 = applyGate(s0, 'X');
  assert.equal(formatKetNotation(stepX0.output).ket, '|1⟩');
  assert.equal(formatExactAmp(stepX0.output.alpha), '0');
  assert.equal(formatExactAmp(stepX0.output.beta), '1');
  assert.equal(getStateAngleDegrees(stepX0.output), 90);

  const stepX1 = applyGate(stepX0.output, 'X');
  assert.equal(formatKetNotation(stepX1.output).ket, '|0⟩');
  assert.equal(formatExactAmp(stepX1.output.alpha), '1');
  assert.equal(formatExactAmp(stepX1.output.beta), '0');
  assert.equal(getStateAngleDegrees(stepX1.output), 0);
});

test('Hadamard maps |0> to |+> at +45 deg and |1> to |-> at -45 deg', () => {
  const h0 = applyGate(createBasisState('0'), 'H');
  assert.equal(formatKetNotation(h0.output).ket, '(1/√2)|0⟩ + (1/√2)|1⟩');
  assert.equal(formatKetNotation(h0.output).alias, '|+⟩');
  assert.equal(formatSignedDecimal(h0.output.alpha), '+0.7071');
  assert.equal(formatSignedDecimal(h0.output.beta), '+0.7071');
  assert.equal(getStateAngleDegrees(h0.output), 45);

  const pH0 = getStateProbabilities(h0.output);
  assert.equal(pH0.p0Fraction, '1/2');
  assert.equal(pH0.p1Fraction, '1/2');
  assert.equal(pH0.p0Percent, '50.0%');
  assert.equal(pH0.p1Percent, '50.0%');
  assert.equal(pH0.sum.toFixed(4), '1.0000');

  const h1 = applyGate(createBasisState('1'), 'H');
  assert.equal(formatKetNotation(h1.output).ket, '(1/√2)|0⟩ - (1/√2)|1⟩');
  assert.equal(formatKetNotation(h1.output).alias, '|-⟩');
  assert.equal(formatSignedDecimal(h1.output.alpha), '+0.7071');
  assert.equal(formatSignedDecimal(h1.output.beta), '-0.7071');
  assert.equal(getStateAngleDegrees(h1.output), -45);

  const pH1 = getStateProbabilities(h1.output);
  assert.equal(pH1.p0, 0.5);
  assert.equal(pH1.p1, 0.5);
  assert.equal(pH1.sum.toFixed(4), '1.0000');
});

test('H(H|0>) = |0> and H(H|1>) = |1> with zero floating-point drift', () => {
  const hh0 = applyGate(applyGate(createBasisState('0'), 'H').output, 'H');
  assert.equal(formatKetNotation(hh0.output).ket, '|0⟩');
  assert.equal(formatExactAmp(hh0.output.alpha), '1');
  assert.equal(formatExactAmp(hh0.output.beta), '0');
  assert.deepEqual(hh0.row0Simplified, ['1/2', '1/2']);
  assert.deepEqual(hh0.row1Simplified, ['1/2', '-1/2']);

  const hh1 = applyGate(applyGate(createBasisState('1'), 'H').output, 'H');
  assert.equal(formatKetNotation(hh1.output).ket, '|1⟩');
  assert.equal(formatExactAmp(hh1.output.alpha), '0');
  assert.equal(formatExactAmp(hh1.output.beta), '1');
});

test('evaluateCircuit tracks per-step states, supports removing middle gates, and toggles initial basis', () => {
  const evalHH = evaluateCircuit('0', ['H', 'H']);
  assert.equal(evalHH.steps.length, 2);
  assert.equal(formatKetNotation(evalHH.steps[0].input).ket, '|0⟩');
  assert.equal(formatKetNotation(evalHH.steps[0].output).alias, '|+⟩');
  assert.equal(formatKetNotation(evalHH.steps[1].input).alias, '|+⟩');
  assert.equal(formatKetNotation(evalHH.steps[1].output).ket, '|0⟩');
  assert.equal(formatKetNotation(evalHH.finalState).ket, '|0⟩');

  // Toggle initial basis from |0> to |1> with [H] on wire
  const evalH1 = evaluateCircuit('1', ['H']);
  assert.equal(formatKetNotation(evalH1.finalState).alias, '|-⟩');
  assert.equal(getStateAngleDegrees(evalH1.finalState), -45);

  // Removing a gate from [H, X, H] -> [H, H] restores |0>
  const withMiddle = evaluateCircuit('0', ['H', 'X', 'H']);
  assert.equal(formatKetNotation(withMiddle.finalState).ket, '|0⟩');
  const afterRemoveFirstH = evaluateCircuit('0', ['X', 'H']);
  assert.equal(formatKetNotation(afterRemoveFirstH.finalState).alias, '|-⟩');
});

test('Identity gate leaves basis and superposition states unchanged', () => {
  const h0 = applyGate(createBasisState('0'), 'H').output;
  const ih0 = applyGate(h0, 'I');
  assert.equal(formatKetNotation(ih0.output).ket, formatKetNotation(h0).ket);
  assert.equal(formatSignedDecimal(ih0.output.alpha), '+0.7071');
  assert.equal(formatSignedDecimal(ih0.output.beta), '+0.7071');
});

test('Measurement simulator collapses 100 shots accurately for deterministic and superposition states', () => {
  const m0 = simulateMeasurements(createBasisState('0'), 100);
  assert.equal(m0.count0, 100);
  assert.equal(m0.count1, 0);
  assert.equal(m0.count0 + m0.count1, 100);

  const m1 = simulateMeasurements(createBasisState('1'), 100);
  assert.equal(m1.count0, 0);
  assert.equal(m1.count1, 100);

  // Deterministic RNG alternating 0.25 and 0.75 should yield exact 50 / 50 for |+>
  let tick = 0;
  const detRng = () => (tick++ % 2 === 0 ? 0.25 : 0.75);
  const h0 = applyGate(createBasisState('0'), 'H').output;
  const mPlus = simulateMeasurements(h0, 100, detRng);
  assert.equal(mPlus.count0, 50);
  assert.equal(mPlus.count1, 50);
  assert.equal(mPlus.deltaP0, 0);
  assert.equal(mPlus.deltaP1, 0);
});

test('Self-contained index.html inline vanilla JS kernel matches all quantum identities and circuit evaluation', () => {
  const htmlPath = path.resolve(process.cwd(), 'index.html');
  const html = fs.readFileSync(htmlPath, 'utf-8');
  const scriptMatch = html.match(/<script>([\s\S]*?)<\/script>/);
  assert.ok(scriptMatch, 'index.html should contain an inline <script> block');

  // Verify self-contained: no external <script src=...>
  assert.equal(/<script[^>]+src=/.test(html), false, 'index.html must not depend on external JS scripts');

  const sandbox: { window: Record<string, unknown>; Math: Math } = {
    window: {},
    Math,
  };
  vm.createContext(sandbox);
  vm.runInContext(scriptMatch[1], sandbox);

  const kernel = sandbox.window.__QUBIT_KERNEL__ as {
    createBasisState: (b: '0' | '1') => unknown;
    applyGate: (s: unknown, g: 'X' | 'H' | 'I') => { output: unknown };
    evaluateCircuit: (
      b: '0' | '1',
      seq: ('X' | 'H' | 'I')[]
    ) => { steps: { input: unknown; output: unknown }[]; finalState: unknown };
    formatKetNotation: (s: unknown) => { ket: string; alias: string | null };
    getStateAngleDegrees: (s: unknown) => number;
    getStateProbabilities: (s: unknown) => { p0: number; p1: number; sum: number };
    simulateMeasurements: (s: unknown, n: number, rng?: () => number) => { count0: number; count1: number };
  };

  assert.ok(kernel, 'inline script should expose __QUBIT_KERNEL__');

  // Verify evaluateCircuit on H(H|0>) = |0> directly in the inline vanilla JS
  const circHH = kernel.evaluateCircuit('0', ['H', 'H']);
  assert.equal(circHH.steps.length, 2);
  assert.equal(kernel.getStateAngleDegrees(circHH.steps[0].output), 45);
  assert.equal(kernel.formatKetNotation(circHH.finalState).ket, '|0⟩');
  assert.equal(kernel.getStateAngleDegrees(circHH.finalState), 0);
  assert.equal(kernel.getStateProbabilities(circHH.finalState).sum.toFixed(4), '1.0000');

  // Verify H|1> = |-> at -45 deg
  const circH1 = kernel.evaluateCircuit('1', ['H']);
  assert.equal(kernel.getStateAngleDegrees(circH1.finalState), -45);
  assert.equal(kernel.formatKetNotation(circH1.finalState).ket, '(1/√2)|0⟩ - (1/√2)|1⟩');

  // Verify X|0> = |1>
  const circX0 = kernel.evaluateCircuit('0', ['X']);
  assert.equal(kernel.formatKetNotation(circX0.finalState).ket, '|1⟩');

  // Verify 100 shots sum to 100
  const shots = kernel.simulateMeasurements(circH1.finalState, 100);
  assert.equal(shots.count0 + shots.count1, 100);
});
