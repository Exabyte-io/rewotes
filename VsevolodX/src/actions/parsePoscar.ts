export function parsePOSCAR(poscarContents: string): string {
    const lines = poscarContents.split(/\r?\n/);
    const elementSymbols = lines[5].split(/\s+/).filter((x) => x);
    const elementCounts = lines[6].split(/\s+/).map(Number);
    const cartesianCoordsStart = 8;
  
    let xyzContents = `${elementCounts.reduce((a, b) => a + b)}\nConverted from POSCAR\n`;
  
    let atomIndex = 0;
    for (let i = 0; i < elementCounts.length; i++) {
      for (let j = 0; j < elementCounts[i]; j++) {
        const coordsLine = lines[cartesianCoordsStart + atomIndex].split(/\s+/).filter((x) => x);
        xyzContents += `${elementSymbols[i]} ${coordsLine.join(' ')}\n`;
        atomIndex++;
      }
    }
  
    return xyzContents;
  }
  