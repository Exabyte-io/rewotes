import { FileUploadButton } from "@components"
import { CrystalInput } from "@types"
import { FC, SetStateAction } from "react"

export const Toolbar: FC<{ setInput: (value: SetStateAction<CrystalInput>) => void }> = ({ setInput }) => {
  const handleFile = (input: string) => {
    setInput((prev) => {
      const newInput = { ...prev }
      newInput["crystalBasis"] = input
      return newInput
    })
  }

  return (
    <div className="absolute top-0 z-30 flex h-20 max-h-20 w-full flex-row items-center justify-between gap-5 bg-dark2 px-10 lg:px-20">
      <p className="whitespace-nowrap font-mozart text-3xl uppercase tracking-widest text-light">Materials Designer</p>
      <FileUploadButton handleFile={handleFile} />
    </div>
  )
}
