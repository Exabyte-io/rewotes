"use client"

import { ChangeEventHandler, FC } from "react"

export const SourceEditor: FC<{ handleEditor: (input: string) => void; input: string }> = ({ handleEditor, input }) => {
  const handleChange: ChangeEventHandler<HTMLTextAreaElement> = async (event) => {
    console.log(event.target.value)
    handleEditor(event.target.value)
  }
  return (
    <div className="flex h-full w-full flex-col gap-5 bg-green-200 p-2 lg:p-10">
      <div className="flex flex-col gap-2">
        <p className="whitespace-nowrap text-xl">Source Editor</p>
        <p className="text-sm">{`Try adding 3D vectors, comma separated, to generate additional cubes. One line per cube to set the following properties: Position<x, y, z>`}</p>
      </div>
      {/* TODO: try monaco editor w/ custom language */}
      <textarea
        value={input}
        onChange={handleChange}
        className="flex grow rounded-sm p-2 font-mono text-sm text-black focus:outline-purple-400"
      ></textarea>
    </div>
  )
}
