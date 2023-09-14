import { parseCSV, parseXYZ } from "@/helpers"
import { Explorer, SourceEditor, Toolbar, Visualizer } from "@components"
import { CrystalInput } from "@types"
import { Inter } from "next/font/google"
import { useState } from "react"
import { Vector3 } from "three"

const inter = Inter({ subsets: ["latin"] })

const latticeTypes: Record<string, string> = {
  // TODO: add more lattice types
  cube: `1,-1,1\n1,1,1\n-1,1,1\n-1,-1,1\n1,-1,1\n\n1,-1,-1\n1,1,-1\n-1,1,-1\n-1,-1,-1\n1,-1,-1\n\n1,-1,1\n1,-1,-1\n\n1,1,1\n1,1,-1\n\n-1,1,1\n-1,1,-1\n\n-1,-1,1\n-1,-1,-1\n\n1,-1,1\n1,-1,-1`,
}

const baseLatticeVectors = parseCSV(latticeTypes.cube)

export default function Home() {
  const [input, setInput] = useState<CrystalInput>({
    crystalBasis: "C 0 0 0",
    a: 2,
    b: 2,
    c: 2,
    BC: 90,
    AC: 90,
    AB: 90,
  })

  const pointVectors = parseXYZ(input.crystalBasis)[1] as Vector3[]

  const { a, b, c } = input
  const xCoord = a > 0 ? a / 2 : 0
  const yCoord = b > 0 ? b / 2 : 0
  const zCoord = c > 0 ? c / 2 : 0

  // set the input somehow when this changes

  const min = new Vector3(-xCoord, -yCoord, -zCoord)
  const max = new Vector3(xCoord, yCoord, zCoord)

  const latticeVectors = baseLatticeVectors.map((v) => {
    const exaggeratedV = v.multiply(max) // required for clamp function to work as intended
    return exaggeratedV.clamp(min, max)
  })

  // TODO: is useMemo even optimizing performance here?
  // const vectors = useMemo(() => {
  //   console.log('useMemo')
  //   return vectorInput.map((v) => {
  //     const exaggeratedV = v.multiply(max) // required for clamp function to work as intended
  //     return exaggeratedV.clamp(min, max)
  //   })
  // }, [min, max])

  const handleEditor = (id: keyof CrystalInput, input: string | number) => {
    setInput((prev) => {
      const newInput = { ...prev }
      if (id === "crystalBasis") {
        newInput[id] = input as string
      } else {
        newInput[id] = input as number
      }
      return newInput
    })
  }

  const [hideExplorer, setHideExplorer] = useState<boolean>(false)

  return (
    <>
      <Toolbar setInput={setInput} />
      <div className="flex min-h-screen w-full bg-dark2 pt-20">
        <div className="flex w-full flex-row">
          <Explorer hide={hideExplorer} setHide={setHideExplorer} />
          <div className="smooth-transition-all flex w-full flex-col lg:flex-row">
            <SourceEditor handleEditor={handleEditor} input={input} />
            <Visualizer latticeVectors={latticeVectors} pointVectors={pointVectors} />
          </div>
        </div>
      </div>
    </>
  )
}
