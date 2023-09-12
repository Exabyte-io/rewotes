import { Explorer, SourceEditor, Visualizer } from "@components"
import { Inter } from "next/font/google"
import { useEffect, useState } from "react"

const inter = Inter({ subsets: ["latin"] })

export default function Home() {
  const [input, setInput] = useState<string>("")

  useEffect(() => {
    console.log(input)
  }, [input])

  const handleEditor = async (input: string) => {
    setInput(input)
  }

  return (
    <div className="flex min-h-screen w-full bg-blue-100 pt-20">
      <div className="flex w-full flex-row">
        <Explorer />
        <SourceEditor handleEditor={handleEditor} />
        <Visualizer input={input} />
      </div>
    </div>
  )
}
