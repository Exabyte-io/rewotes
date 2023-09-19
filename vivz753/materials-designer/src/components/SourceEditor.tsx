import { CrystalInput } from "@types"
import clsx from "clsx"
import { ChangeEventHandler, FC, useState } from "react"

const inputIds = ["a", "b", "c", "BC", "AC", "AB"]

export const SourceEditor: FC<{
  handleEditor: (id: keyof CrystalInput, input: string) => void
  input: CrystalInput
}> = ({ handleEditor, input }) => {
  const handleChange: ChangeEventHandler<HTMLTextAreaElement | HTMLInputElement> = async (event) => {
    handleEditor(event.target.id as keyof CrystalInput, event.target.value)
  }

  const [collapseExplorer, setCollapseExplorer] = useState<boolean>(false)

  const toggleEditor = () => {
    setCollapseExplorer((prev) => !prev)
  }

  // orthorhombic, tetragonal, cubic -- can NOT edit angles

  return (
    <div
      className={clsx(
        "relative flex w-full flex-col gap-8 rounded-sm border border-dark1 bg-dark2 px-5 py-8 font-mono text-light focus-within:border focus-within:border-accent lg:px-10 lg:py-10",
        // TODO: fix collpase function
        collapseExplorer ? "h-48" : "h-full"
      )}
    >
      {/* <button
        onClick={toggleEditor}
        className="absolute right-0 top-0 m-4 flex h-6 w-6 items-center justify-center rounded-full bg-accent lg:hidden"
      >
        v
      </button> */}
      <p className="font-mozart text-xl uppercase tracking-widest">Editor</p>
      <div className={clsx("flex w-full flex-col gap-2", collapseExplorer && "hidden")}>
        <p className="whitespace-nowrap text-xl">Crystal Lattice</p>
        {/* TODO: add input dropdown w/ types for crystal structures */}
        <div className="flex flex-wrap">
          {inputIds.map((id) => (
            <div className="my-1 flex w-1/3 flex-col items-start lg:my-2">
              <label className="mx-1 mb-0.5 font-light">{id}</label>
              <input
                disabled={id !== "a" && id !== "b" && id !== "c"}
                id={id}
                value={input[id as keyof CrystalInput]}
                onChange={handleChange}
                className="w-20 rounded-sm p-1.5 text-center font-mono text-sm text-dark2 focus:outline-accent"
              />
            </div>
          ))}
        </div>
      </div>
      <div className={clsx("flex h-full w-full flex-col gap-5", collapseExplorer && "hidden")}>
        <div className="flex flex-col gap-2">
          <p className="whitespace-nowrap text-xl">Crystal Basis</p>
          <p className="text-sm font-light italic">{`Try adding a chemical symbol with a 3D vector to generate additional points.`}</p>
        </div>
        <textarea
          id="crystalBasis"
          placeholder="<symbol x y z> (i.e. C 0 0 0)"
          value={input.crystalBasis}
          onChange={handleChange}
          className="flex grow rounded-sm p-2 font-mono text-sm text-dark2 focus:outline-accent"
        ></textarea>
      </div>
    </div>
  )
}
