import clsx from "clsx"
import { ChangeEventHandler, FC, useState } from "react"

const inputIds = ["a", "b", "c", "BC", "AC", "AB"]

export const SourceEditor: FC<{
  handleEditor: (id: string, input: string) => void
  input: Record<string, string | number>
}> = ({ handleEditor, input }) => {
  const handleChange: ChangeEventHandler<HTMLTextAreaElement | HTMLInputElement> = async (event) => {
    handleEditor(event.target.id, event.target.value)
  }

  const [collapseExplorer, setCollapseExplorer] = useState<boolean>(false)

  const toggleEditor = () => {
    setCollapseExplorer((prev) => !prev)
  }

  return (
    <div
      className={clsx(
        "relative flex w-full flex-col gap-10 bg-green-200 px-5 py-10 lg:p-10",
        collapseExplorer ? "h-48" : "h-full"
      )}
    >
      <button
        onClick={toggleEditor}
        className="absolute right-0 top-0 m-4 flex h-6 w-6 items-center justify-center rounded-full bg-purple-500 lg:hidden"
      >
        v
      </button>
      <div className={clsx("flex h-full w-full flex-col gap-5", collapseExplorer && "hidden")}>
        <div className="flex flex-col gap-2">
          <p className="whitespace-nowrap text-xl">Crystal Lattice</p>
        </div>
        {/* TODO: try monaco editor w/ custom language */}
        <div className="flex flex-wrap">
          {inputIds.map((id) => (
            <div className="flex w-1/3  flex-col items-center">
              <label>{id}</label>
              <input
                id={id}
                value={input[id]}
                onChange={handleChange}
                className="w-20 rounded-sm p-2 font-mono text-sm text-black focus:outline-purple-400"
              />
            </div>
          ))}
        </div>
      </div>
      <div className={clsx("flex h-full w-full flex-col gap-5", collapseExplorer && "hidden")}>
        <div className="flex flex-col gap-2">
          <p className="whitespace-nowrap text-xl">Crystal Basis</p>
          <p className="text-sm">{`Try adding 3D vectors, comma separated, to generate additional cubes. One line per cube to set the following properties: Position<x, y, z>`}</p>
        </div>
        {/* TODO: try monaco editor w/ custom language */}
        <textarea
          id="crystalBasis"
          value={input.crystalBasis}
          onChange={handleChange}
          className="flex grow rounded-sm p-2 font-mono text-sm text-black focus:outline-purple-400"
        ></textarea>
      </div>
    </div>
  )
}
