import React, { useRef, useState } from 'react';

const Sphere = (props) => {
  // This reference will give us direct access to the mesh
  const mesh = useRef();

  // Set up state for the hovered and active state
  const [hovered, setHover] = useState(false);

  return (
    <mesh
      {...props}
      ref={mesh}
      scale={[0.25, 0.25, 0.25]}
      onPointerOver={(e) => {
        props.setHoveredElement(props.id);
        setHover(true);
      }}
      onPointerOut={(e) => {
        props.setHoveredElement(null);
        setHover(false);
      }}
    >
      <sphereGeometry args={[1, 64, 64]} />
      <meshStandardMaterial color={hovered ? 'skyblue' : 'green'} />
    </mesh>
  );
};

export default Sphere;
