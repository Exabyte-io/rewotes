import { extend, useThree } from 'react-three-fiber';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';
import React from 'react';

extend({ OrbitControls });

const Controls = () => {
  const { camera, gl } = useThree();

  return (
    <orbitControls

      args={[camera, gl.domElement]}
    />
  );
};

export default Controls;
