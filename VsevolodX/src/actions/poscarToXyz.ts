export function poscarToXYZ(contents: string): string | null {

  const lines = contents.split('\n');

  const comment = lines[0];
  const scaleStr = lines[1].trim().split(' ');
  const scale = parseFloat(scaleStr[0]);
  const latticeVectors = lines.slice(2, 5).map((line) => line.split(' ').map(parseFloat));
  const elementSymbols = lines[5].trim().split(' ');
  const elementCounts = lines[6].trim().split(' ').map((x) => parseInt(x, 10));
  const cooridinatesType = lines[7][0]; //first character designates Cartesian or Direct coordinates 
  const positions = lines.slice(8, 8 + sum(elementCounts)).map((line) => line.split(' ').map(parseFloat));

  let elements: string[] = [];
  elementSymbols.forEach((symbol, idx) => {
    for (let i = 0; i < elementCounts[idx]; i++) {
      elements.push(symbol);
    }
  });
  
  const convertedPositions = positions.map((pos) => {
    const [x, y, z] = pos;
    const latticeX = [
      x * latticeVectors[0][0],
      x * latticeVectors[1][0],
      x * latticeVectors[2][0],
    ];
    const latticeY = [
      y * latticeVectors[0][1],
      y * latticeVectors[1][1],
      y * latticeVectors[2][1],
    ];
    const latticeZ = [
      z * latticeVectors[0][2],
      z * latticeVectors[1][2],
      z * latticeVectors[2][2],
    ];
    return [
      scale * (latticeX[0] + latticeY[0] + latticeZ[0]),
      scale * (latticeX[1] + latticeY[1] + latticeZ[1]),
      scale * (latticeX[2] + latticeY[2] + latticeZ[2]),
    ];
  });

  const data = { elements, convertedPositions };
  const numAtoms = data.elements.length;
  let output = `${numAtoms}\n`
  output += `${comment}. Converted from POSCAR\n`;

  data.elements.forEach((element, idx) => {
    const position = data.convertedPositions[idx].join(' ');
    output += `${element} ${position}\n`;
  });

    return output;
  }

function sum(arr: number[]): number {
    return arr.reduce((a, b) => a + b, 0);
  }