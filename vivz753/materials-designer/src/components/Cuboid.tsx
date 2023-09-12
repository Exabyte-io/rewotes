import React, { useRef } from "react"
import { useLoader, useFrame, extend } from "@react-three/fiber"
import { MeshWobbleMaterial } from "@react-three/drei"
extend({MeshWobbleMaterial})
import { Mesh } from "three"

export const Cuboid = () => {
  const mesh = useRef<Mesh>(null!)
  useFrame(() => (mesh.current.rotation.x = mesh.current.rotation.y += 0.01))

  return (
    <mesh castShadow ref={mesh}>
      <boxBufferGeometry attach="geometry" args={[1, 2, 1]} />
      <MeshWobbleMaterial attach="material" color="red" />
    </mesh>
  )
}
