import React, { useRef, useState, useEffect } from 'react';
import * as THREE from 'three';

const Line = ({ connection }) => {
  const [points, setPoints] = useState([]);
  console.log(connection)
  // This reference will give us direct access to the mesh
  const mesh = useRef();
  useEffect(() => {
    const vectors = [];
    connection.forEach((element) => {
      console.log(element);
      vectors.push(new THREE.Vector3(element.x, element.y, element.z));
    });
    setPoints([...vectors]);
  }, []);

  const lineGeometry = new THREE.BufferGeometry().setFromPoints(points);

  return (
    <line ref={mesh} geometry={lineGeometry}>
      <lineBasicMaterial
        attach='material'
        color={'#9c88ff'}
        linewidth={100}
        linecap={'round'}
        linejoin={'round'}
      />
    </line>
  );
};

export default Line;
