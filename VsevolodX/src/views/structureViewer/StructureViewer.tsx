import { Card } from '@blueprintjs/core'
import React, { useEffect, useState } from 'react'
import styles from './StructureViewer.module.scss';
import { Canvas } from '@react-three/fiber';
import { Sphere, OrbitControls } from '@react-three/drei';
import { Object3D, Vector3 } from 'three';
import SourceContext from '../../context/SourceContext';
import ViewHeading from '../../components/view_heading/ViewHeading';

import { useSettings } from '../../context/SettingsContext';
import VStack from '../../components/utils/VStack';


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
  const atomsData = settings.settings.atomsData || [];
  console.log(atomsData); //TODO: REMOVE before production

  const viewBoxSize: number = 20;
  const [atoms, setAtoms] = useState<Atom[]>([]); //TODO: save atoms to a Context

  useEffect(() => {
    if (isValidXYZFormat) {
      if (!parseAndRender(source, setAtoms)) {
        setIsValidXYZFormat(false);
      }
    }
  }, [source, isValidXYZFormat, setIsValidXYZFormat]);

  Object3D.DEFAULT_UP.set(0, 0, 1);

  return (
    <Card className={styles.StructureViewer}>
      <VStack>
      <ViewHeading>
        <h4>Structure Viewer</h4>
      </ViewHeading>
    <Canvas className={styles.Canvas + ' ' + theme}>
      <OrbitControls />
      <ambientLight />
      <pointLight position={[viewBoxSize, viewBoxSize, viewBoxSize]} />
      <gridHelper args={[viewBoxSize, viewBoxSize, '#CCCCCC', '#AAAAAA']} rotation-x={Math.PI / 2} />
      <axesHelper args={[viewBoxSize/2]} />
      {atoms.map((atom, index) => {
        const colorObj = atomsData.find((color) => color.element === atom.element);
        const defaultColor = settings.settings.defaultAtomColor;
        const materialColor = colorObj ? colorObj.color : defaultColor;
        return (
          <Sphere key={atom.position.x} args={[0.1]} position={atom.position.toArray()}>
          <meshStandardMaterial attach="material" color={materialColor} />
        </Sphere>
      )}
      )};
  
    </Canvas>
      </VStack>
    </Card>
  );
};

export default StructureViewer