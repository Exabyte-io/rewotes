import { Card } from '@blueprintjs/core'
import React, { useEffect, useState } from 'react'
import styles from './StructureViewer.module.scss';
import { Canvas } from '@react-three/fiber';
import { Sphere, OrbitControls } from '@react-three/drei'
import { Vector3 } from 'three';
import SourceContext from './SourceContext';


interface Atom {
  element: string;
  position: Vector3;
}
const parseXYZ = (text: string): Atom[] => {
  const lines = text.split('\n');
  const atomCount = parseInt(lines[0].trim(), 10);
  const atoms: Atom[] = [];

  for (let i = 2; i < atomCount + 2; i++) {
    const line = lines[i].split(/\s+/);
    const element = line[0];
    const x = parseFloat(line[1]);
    const y = parseFloat(line[2]);
    const z = parseFloat(line[3]);
    atoms.push({ element, position: new Vector3(x, y, z) });
  }

  return atoms;
};

function parseAndRender(input: string, setAtoms: { (value: React.SetStateAction<Atom[]>): void; (arg0: Atom[]): void; }) {
  try {

    const parsedAtoms = parseXYZ(input);
    setAtoms(parsedAtoms); 
    console.log('Pasing successfull');
    return true;
  } catch (error) {
    console.error('Failed at parsing input for XYZ:', error);
    return false;
  }
}

const StructureViewer: React.FC = () => {
  const { source, isValidXYZFormat, setIsValidXYZFormat } = React.useContext(SourceContext);
  
  const [atoms, setAtoms] = useState<Atom[]>([]);

  useEffect(() => {
    if (isValidXYZFormat) {
      if (!parseAndRender(source, setAtoms)) {
        setIsValidXYZFormat(false);
      }
    }
  }, [source, isValidXYZFormat, setIsValidXYZFormat]);


  return (
    <Card className={styles.StructureViewer}>
    <Canvas style={{ width: '100%', height: '100%' }}>
      <OrbitControls />
      <ambientLight />
      <pointLight position={[10, 10, 10]} />
      {atoms.map((atom, index) => (
        <Sphere key={index} args={[0.1]} position={atom.position.toArray()}>
          <meshStandardMaterial attach="material" color={atom.element === 'Si' ? 'red' : 'blue'} />
        </Sphere>
      ))}
    </Canvas>
    </Card>
  );
};

export default StructureViewer