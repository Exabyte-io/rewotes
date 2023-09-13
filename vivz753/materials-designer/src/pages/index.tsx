import { Explorer, SourceEditor, Visualizer } from "@components"
import { Inter } from "next/font/google"
import { useEffect, useState } from "react"

const inter = Inter({ subsets: ["latin"] })

export default function Home() {
  const [input, setInput] = useState<string>(`1,0,0\n0,1,0\n-1,0,0\n0,-1,0\n1,0,0`)

  useEffect(() => {
    console.log(input)
  }, [input])

  const handleEditor = async (input: string) => {
    setInput(input)
  }

  const [hideExplorer, setHideExplorer] = useState<boolean>(false)

  return (
    <div className="flex min-h-screen w-full bg-blue-100 pt-20">
      <div className="flex w-full flex-row">
        <Explorer hide={hideExplorer} setHide={setHideExplorer} />
        <div className="smooth-transition-all flex w-full flex-col lg:flex-row">
          <SourceEditor handleEditor={handleEditor} input={input} />
          <Visualizer input={input} />
        </div>
      </div>
    </div>
  )
}
