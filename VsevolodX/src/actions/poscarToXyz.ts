export function poscarToXYZ(contents: string): string | null {
  
  const lines = contents.split('\n');

  console.log('past regex test')
  const scale = parseFloat(lines[1].trim());
  const latticeVectors = lines.slice(2, 5).map((line) => line.split(' ').map(parseFloat));
  const elementSymbols = lines[5].trim().split(' ');
  const elementCounts = lines[6].trim().split(' ').map((x) => parseInt(x, 10));
  const positions = lines.slice(8, 8 + sum(elementCounts)).map((line) => line.split(' ').map(parseFloat));

  let elements: string[] = [];
  elementSymbols.forEach((symbol, idx) => {
    for (let i = 0; i < elementCounts[idx]; i++) {
      elements.push(symbol);
    }
  });
  
  const data = { elements, positions };
  const numAtoms = data.elements.length;
  let output = `${numAtoms}\nConverted from POSCAR\n`;

  data.elements.forEach((element, idx) => {
    const position = data.positions[idx].join(' ');
    output += `${element} ${position}\n`;
  });

    return output;
  }

function sum(arr: number[]): number {
    return arr.reduce((a, b) => a + b, 0);
  }