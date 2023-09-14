import { Explorer, SourceEditor, Visualizer } from "@components"
import { CrystalInput } from "@types"
import { Inter } from "next/font/google"
import { useState } from "react"
import { Vector3 } from "three"

const inter = Inter({ subsets: ["latin"] })

const crystalTypes: Record<string, string> = {
  cube: `1,-1,1\n1,1,1\n-1,1,1\n-1,-1,1\n1,-1,1\n\n1,-1,-1\n1,1,-1\n-1,1,-1\n-1,-1,-1\n1,-1,-1\n\n1,-1,1\n1,-1,-1\n\n1,1,1\n1,1,-1\n\n-1,1,1\n-1,1,-1\n\n-1,-1,1\n-1,-1,-1\n\n1,-1,1\n1,-1,-1`,
}

const parseInput = (input: string): Vector3[] => {
  const lines = input.split("\n")
  const tokens = lines.map((s) => s.split(","))
  const cleanedTokens = tokens.map((line) => new Vector3(parseFloat(line[0]), parseFloat(line[1]), parseFloat(line[2])))
  return cleanedTokens
}

export default function Home() {
  const [input, setInput] = useState<CrystalInput>({
    crystalBasis: crystalTypes["cube"],
    a: 2,
    b: 2,
    c: 2,
    BC: 90,
    AC: 90,
    AB: 90,
  })
  console.log("redner")

  const vectorInput = parseInput(input.crystalBasis)

  const { a, b, c } = input
  const xCoord = a > 0 ? a / 2 : 0
  const yCoord = b > 0 ? b / 2 : 0
  const zCoord = c > 0 ? c / 2 : 0

  // set the input somehow when this changes

  const min = new Vector3(-xCoord, -yCoord, -zCoord)
  const max = new Vector3(xCoord, yCoord, zCoord)

  const vectors = vectorInput.map((v) => {
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

  const handleEditor = async (id: keyof CrystalInput, input: string | number) => {
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
    <div className="flex min-h-screen w-full bg-blue-100 pt-20">
      <div className="flex w-full flex-row">
        <Explorer hide={hideExplorer} setHide={setHideExplorer} />
        <div className="smooth-transition-all flex w-full flex-col lg:flex-row">
          <SourceEditor handleEditor={handleEditor} input={input} />
          <Visualizer vectors={vectors} />
        </div>
      </div>
    </div>
  )
}
