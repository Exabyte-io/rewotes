import React, { useState } from 'react';
import { Canvas } from 'react-three-fiber';
import Sphere from './Sphere';
import Controls from './Controls';
import Line from './Connect';
import styled from 'styled-components';

const CanvasArea = ({ elements, connections, setHoveredElement }) => {
  console.log('Connections', connections);
  return (
    <CanvasStyles>
      <Canvas camera={{ position: [5, 5, 2] }}>
        <gridHelper />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {elements.map((element) => (
          <Sphere
            id={element.id}
            key={element.id}
            position={[element.x, element.y, element.z]}
            setHoveredElement={setHoveredElement}
          />
        ))}
        {connections.map((connection, idx) => (
          <Line key={idx} connection={connection} />
        ))}
        <Controls />
      </Canvas>
    </CanvasStyles>
  );
};

const CanvasStyles = styled.div`
  height: 100vh;
  width: 70%;
`;

export default CanvasArea;
