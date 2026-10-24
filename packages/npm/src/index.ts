export { MERSENNE_61, MERSENNE_127, mersenneMod, mersenneMul, phi, PolynomialHash, HashProfile } from './polynomial.js';
export type { Integer } from './polynomial.js';
export { Leaf, Internal, RepeatNode, ropeLen, ropeHash, ropeHeight, ropeFromBytes, ropeToBytes, validateRope, structuralEquals, RopeValidationError } from './nodes.js';
export type { Node, NonEmptyNode } from './nodes.js';
export { ALPHA_NUM, ALPHA_DEN, ropeConcat, ropeSplit, ropeRepeat, ropeSubstrHash, ropeInsert, ropeDelete, contentEquals } from './rope.js';
export { SlidingWindow } from './sliding.js';
/** npm implementation version. */
export const VERSION = '0.3.3';
/** The independently pinned upstream behavior release. */
export const PYTHON_VERSION = '0.2.2';
export { ropeEqual, ropeLcp, ropeRangeEqual, ropeRangeLcp } from './comparison.js';
