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
      <sphereGeometry args={[1, 32, 32]} />
      <meshStandardMaterial color={hovered ? '#2471a3' : '#229954'} />
    </mesh>
  );
};

export default Sphere;
