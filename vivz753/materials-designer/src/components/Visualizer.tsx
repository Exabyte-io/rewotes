"use client"
import { Line } from "@components"
import { TrackballControls } from "@react-three/drei"
import { Canvas } from "@react-three/fiber"
import { FC } from "react"
import { Vector3 } from "three"
import Box from "./Box"

const parseInput = (input: string): Vector3[] => {
  const lines = input.split("\n")
  const tokens = lines.map((s) => s.split(","))
  const cleanedTokens = tokens.map((line) => new Vector3(parseFloat(line[0]), parseFloat(line[1]), parseFloat(line[2])))
  console.log("lines", lines)
  console.log("tokens", tokens)
  console.log("cleanedTokens", cleanedTokens)
  return cleanedTokens
}

export const Visualizer: FC<{ input: string }> = ({ input }) => {
  const vectors = parseInput(input)

  return (
    <div className="flex h-full w-full grow flex-col justify-center">
      <p className="text-center text-xl">Visualizer</p>
      <Canvas className="h-full w-full bg-pink-200">
        <TrackballControls />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {vectors.map((vector, i) => (
          <Box key={i} position={vector} />
        ))}
        <Line vectors={vectors} />
        {/* <Box position={[-1.2, 0, 0]} /> */}
        {/* <Box position={[1.2, 0, 0]} /> */}
      </Canvas>
    </div>
  )
}
