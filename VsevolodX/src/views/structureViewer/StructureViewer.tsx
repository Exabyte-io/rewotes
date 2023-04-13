import { Card } from '@blueprintjs/core'
import React, { useEffect, useState } from 'react'
import styles from './StructureViewer.module.scss';
import { Canvas } from '@react-three/fiber';
import { Sphere, OrbitControls } from '@react-three/drei'
import { Vector3 } from 'three';
import SourceContext from '../../context/SourceContext';
import ViewHeading from '../../components/view_heading/ViewHeading';

import settingsJSON from './temp.json'; //TODO: remoove this temp color setting and move to global context
import { useSettings } from '../../context/SettingsContext';

const atomColors = settingsJSON.atomColors;
interface Atom {
  element: string;
  position: Vector3;
}

//TODO: add solution for the case when first 2 lines are empty 
const parseXYZ = (text: string): Atom[] => {
  const lines: string[] = text.split('\n');
  const atomCount: number = parseInt(lines[0].trim(), 10); // first value in XYZ file is number of atoms
  const atoms: Atom[] = [];
  const comment: string = lines[1]; // second is comment

  for (let i = 2; i < atomCount + 2; i++) {
    const line = lines[i].split(/\s+/); // splitting "element x y z" into "element", "x", "y", "z"
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
    console.log('Pasing successfull'); //TODO: REMOVE before production
    return true;
  } catch (error) {
    console.error('Failed at parsing input for XYZ:', error);
    return false;
  }
}

const StructureViewer: React.FC = () => {
  const settings = useSettings();
  const theme = settings.settings.theme;
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
      <ViewHeading>
        <h4>Structure Viewer</h4>
      </ViewHeading>
    <Canvas className={styles.Canvas}>
      <OrbitControls />
      <ambientLight />
      <pointLight position={[10, 10, 10]} />
      {atoms.map((atom, index) => {
         const colorObj = atomColors.find((color) => color.element === atom.element);
         const materialColor = colorObj ? colorObj.color : 'white';
      return (
        <Sphere key={index} args={[0.1]} position={atom.position.toArray()}>
          <meshStandardMaterial attach="material" color={materialColor} />
        </Sphere>
      )}
      )};

    </Canvas>
    </Card>
  );
};

export default StructureViewer