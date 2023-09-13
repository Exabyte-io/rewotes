import { FC, useRef } from "react"
import { Mesh } from "three"

export const Sphere: FC = () => {
  const meshRef = useRef<Mesh>(null!)
  return (
    <mesh ref={meshRef}>
      <sphereGeometry args={[2, 32, 16]} />
      {/* <meshStandardMaterial color="blue" /> */}
      <meshBasicMaterial color="orange" />
    </mesh>
  )
}
