import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import ts from 'typescript';
const upstream = process.argv[2];
if (!upstream) throw new Error('Usage: node scripts/audit-retained-arithmetic.mjs <authoritative-0.2.2-directory>');
const printer = ts.createPrinter({removeComments: true});
const digest = s => createHash('sha256').update(s).digest('hex');
function collect(text) {
  const file = ts.createSourceFile('polynomial.ts', text, ts.ScriptTarget.ES2022, true);
  const result = new Map();
  for (const statement of file.statements) {
    if (ts.isFunctionDeclaration(statement)) result.set(statement.name.text, printer.printNode(ts.EmitHint.Unspecified, statement, file));
    if (ts.isClassDeclaration(statement) && statement.name.text === 'PolynomialHash') {
      for (const member of statement.members) if (ts.isMethodDeclaration(member)) result.set(`PolynomialHash.${member.name.getText(file)}`, printer.printNode(ts.EmitHint.Unspecified, member, file));
    }
  }
  return result;
}
const old = collect(readFileSync(`${upstream}/src/polynomial.ts`, 'utf8'));
const current = collect(readFileSync(new URL('../src/polynomial.ts', import.meta.url), 'utf8'));
const rows = [...old].map(([name, source]) => ({name, ast_sha256: digest(source), unchanged: source === current.get(name)}));
if (!rows.every(row => row.unchanged)) throw new Error('Retained arithmetic implementation changed');
const result = {baseline: 'published npm 0.2.2 authoritative source', methodology: 'TypeScript printer-normalized function/method declarations, comments removed; constructor/profile additions excluded and covered separately by regressions', functions: rows, passed: true};
writeFileSync(new URL('../evidence/retained-arithmetic-audit.json', import.meta.url), JSON.stringify(result, null, 2) + '\n');
console.log(`All ${rows.length} existing polynomial functions/methods have identical normalized source`);
