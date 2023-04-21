import { poscarToXYZ } from './poscarToXyz';

describe('poscarToXYZ', () => {
  it('should correctly convert POSCAR to XYZ format', () => {
    const poscarInput = `POSCAR test
    1.0
1.00000000000000 0.00000000000000 0.00000000000000
0.00000000000000 2.00000000000000 0.00000000000000
0.00000000000000 0.00000000000000 3.00000000000000
H O
2 1
Direct
1.00000000000000 0.00000000000000 0.00000000000000
1.00000000000000 0.50000000000000 0.00000000000000
-1.00000000000000 0.50000000000000 1.00000000000000
`;

    const expectedOutput = `3
POSCAR test. Converted from POSCAR
H 1 0 0
H 1 1 0
O -1 1 3
`;

    const result = poscarToXYZ(poscarInput);
    expect(result).toEqual(expectedOutput);
  });

  it('should return error for invalid input', () => {
    const invalidInput = `Invalid POSCAR content`;

    const result = poscarToXYZ(invalidInput);
    expect(result).toHaveErrorMessage('Invalid input');
  });
});
