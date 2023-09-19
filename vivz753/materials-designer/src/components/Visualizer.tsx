"use client"
import { Line, Sphere } from "@components"
import { TrackballControls } from "@react-three/drei"
import { Canvas } from "@react-three/fiber"
import { FC } from "react"
import { Vector3 } from "three"

const colors = ["orange"]

export const Visualizer: FC<{ latticeVectors: Vector3[][]; pointVectors: Vector3[] }> = ({
  latticeVectors,
  pointVectors,
}) => {
  return (
    <div className="relative flex h-full w-full flex-col justify-center rounded-sm border border-dark1">
      <Canvas className="h-full w-full bg-dark2">
        <TrackballControls />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {pointVectors.map((vector, i) => {
          const index = i % colors.length
          return <Sphere key={i} position={vector} color={colors[index]} />
        })}
        {latticeVectors.map((line) => (
          <Line vectors={line} />
        ))}
        <Line color="red" vectors={[new Vector3(0, 0, 0), new Vector3(0.5, 0, 0)]} />
        <Line color="yellow" vectors={[new Vector3(0, 0, 0), new Vector3(0, 0.5, 0)]} />
        <Line color="purple" vectors={[new Vector3(0, 0, 0), new Vector3(0, 0, 0.5)]} />
      </Canvas>
      <p className="absolute top-0 mt-10 px-10 text-center font-mozart text-xl uppercase tracking-widest text-light">
        Visualizer
      </p>
      <div className="absolute bottom-0 lg:bottom-auto right-0 m-2 border border-dark1 p-2 font-mozart text-xl uppercase tracking-widest text-light lg:top-0 lg:m-10 lg:p-5">
        <div className="flex flex-col">
          <p>Legend</p>
          <div className="flex flex-row items-center gap-2">
            <span className="h-2 w-2 bg-red-500" />
            <p>x</p>
          </div>
          <div className="flex flex-row items-center gap-2">
            <span className="h-2 w-2 bg-yellow-500" />
            <p>y</p>
          </div>
          <div className="flex flex-row items-center gap-2">
            <span className="h-2 w-2 bg-purple-500" />
            <p>z</p>
          </div>
          <div className="flex flex-row items-center gap-2">
            <span className="h-[1px] w-2 bg-accent" />
            <p>lattice</p>
          </div>
          <div className="flex flex-row items-center gap-2">
            <span className="h-2 w-2 rounded-full bg-orange-500" />
            <p>atom</p>
          </div>
        </div>
      </div>
    </div>
  )
}
