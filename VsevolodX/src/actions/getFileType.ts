type FileType = 'xyz' | 'poscar' | 'unknown';

export function getFileType(contents: string): FileType{
  const lines = contents.split(/\r?\n/);
  if (!lines.length) {
    return 'unknown';
  }

  const firstLine = lines[0].trim();
  const firstLineAsNumber = Number(firstLine);

  if (!isNaN(firstLineAsNumber) && lines.length > 2 && lines[2].split(/\s+/).length >= 4) {
    return 'xyz';
  }

  if (lines.length > 6 && (lines[6].toLowerCase().startsWith("d") || lines[6].toLowerCase().startsWith("c"))) {
    return 'poscar';
  }

  return 'unknown';
}


export function isPOSCAR(contents: string): boolean {
    const lines = contents.split('\n');
    if (lines.length < 8) return false;
    
    const scaleFactor = parseFloat(lines[1]);
    if (isNaN(scaleFactor)) return false;
    
    for (let i = 2; i < 5; i++) {
      const vectorComponents = lines[i].split(/\s+/).filter((s) => s);
      if (vectorComponents.length !== 3) return false;
      for (const component of vectorComponents) {
        if (isNaN(parseFloat(component))) return false;
      }
    }
    
    const elements = lines[5].split(/\s+/).filter((s) => s);
    if (elements.length === 0) return false;
    
    const counts = lines[6].split(/\s+/).map(Number);
    if (counts.length !== elements.length) return false;
  
    return true;
  }
  
function isXYZ(contents: string): boolean {
  const lines = contents.split('\n');
  if (lines.length < 3) return false;

  const atomCount = parseInt(lines[0].trim());
  if (isNaN(atomCount) || atomCount <= 0) return false;

  return /^\s*\d+\s+.+/.test(lines[2]);
}
