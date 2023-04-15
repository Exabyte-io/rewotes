import { Card } from '@blueprintjs/core'
import React, { useEffect, useState } from 'react'
import styles from './StructureViewer.module.scss';
import { Canvas } from '@react-three/fiber';
import { Sphere, OrbitControls } from '@react-three/drei';
import { Object3D, Vector3 } from 'three';
import ViewHeading from '../../components/view_heading/ViewHeading';

import { useSettings } from '../../context/SettingsContext';
import VStack from '../../components/utils/VStack';
import { useAtomsContext } from '../../context/AtomsContext';

const StructureViewer: React.FC = () => {
  const settings = useSettings();
  const theme = settings.settings.theme;
  const atomsDisplayData = settings.settings.atomsDisplayData || [];

  const viewBoxSize: number = 20;
  const {atoms, updateAtoms} = useAtomsContext(); //TODO: save atoms to a Context
  console.log('viewr atoms: ,', atoms);
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
        const colorObj = atomsDisplayData.find((color) => color.element === atom.element);
        const defaultColor = settings.settings.defaultAtomColor;
        const materialColor = colorObj ? colorObj.color : defaultColor;
        return (
          <Sphere key={atom.id} args={[0.1]} position={atom.position.toArray()}>
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