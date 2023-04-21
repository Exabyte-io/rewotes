import { Vector3 } from "three";
import { Atom } from "../context/AtomsContext";

type ParsedXYZResult = {
    isValid: boolean;
    atoms?: Atom[];
  };

export default function parseXYZ (text: string): ParsedXYZResult {
    const FIRST_LINES = 2;
    const lines = text.split('\n');
    const atomCount = parseInt(lines[0].trim(), 10);
    const newAtoms: Atom[] = [];
  
    if (isNaN(atomCount) || lines.length < atomCount + FIRST_LINES) {
      return { isValid: false }; //Must be the exact amount of lines
    }
  
    for (let i = FIRST_LINES; i < atomCount + FIRST_LINES; i++) {
      const line = lines[i].split(/\s+/);
      if (line.length < 4) { // "element, x, y, z" -- exactly 4 objects
        return { isValid: false };
      }
      const element = line[0]; //TODO: Add validation against elements enum
      
      for (let j = 1; j < line.length; j++) {
        const elementValue = line[j];
        if (!/^-?\d*\.?\d+$/.test(elementValue)) {
          return { isValid: false };
        }
      }
      const x = parseFloat(line[1]);
      const y = parseFloat(line[2]);
      const z = parseFloat(line[3]);
      if (isNaN(x) || isNaN(y) || isNaN(z)) {
        return { isValid: false };
      }
  
      newAtoms.push({
        id: i - FIRST_LINES,
        element: element,
        position: new Vector3(x, y, z),
      });
    }
  
    return { isValid: true, atoms: newAtoms };
  };