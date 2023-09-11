'use client'

import { Canvas } from "@react-three/fiber"
import { useMemo } from 'react';

export default function Home() {
  return (
    <div id="canvas-container">
      <Canvas>
        {/* <ambientLight intensity={0.1} />
        <directionalLight color="red" position={[0, 0, 5]} />
        <mesh>
          <boxGeometry />
          <meshStandardMaterial />
        </mesh> */}
      </Canvas>
    </div>
  )
}
