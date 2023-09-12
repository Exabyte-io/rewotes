"use client"
import { Canvas } from "@react-three/fiber"
import { useRef, FC } from "react"
import { Mesh } from "three"
import { OrbitControls } from "@react-three/drei"
// import { Cuboid } from "@components"
import Box from "./Box"

const parseInput = (input: string): GLfloat[][] => {
  const lines = input.split('\n')
  const tokens = lines.map((s) => s.split(','))
  const cleanedTokens = tokens.map((line) => line.map((token) => parseFloat(token.trim())))
  console.log('lines', lines)
  console.log('tokens', tokens)
  console.log('cleanedTokens', cleanedTokens)
  return cleanedTokens
}

export const Visualizer: FC<{ input: string }> = ({ input }) => {
  const vectors = parseInput(input)

  return (
    <div className="flex h-full w-full grow flex-col justify-center">
      <p className="text-center text-xl">Visualizer</p>
      <Canvas className="h-full w-full bg-pink-200">
        <OrbitControls enableZoom={false} />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {vectors.map((vector, i) => (
          <Box key={i} position={vector} />
        ))}
        <Box position={[-1.2, 0, 0]} />
        <Box position={[1.2, 0, 0]} />
        {/* <Cuboid /> */}
        {/* <Sphere /> */}
      </Canvas>
    </div>
  )
}

const Sphere = () => {
  const meshRef = useRef<Mesh>(null!)
  return (
    <mesh ref={meshRef}>
      <sphereGeometry args={[2, 32, 16]} />
      {/* <meshStandardMaterial color="blue" /> */}
      <meshBasicMaterial color="orange" />
    </mesh>
  )
}
