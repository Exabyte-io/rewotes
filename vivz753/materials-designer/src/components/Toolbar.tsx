import { FileUploadButton } from "@components"
import { FC } from "react"

export const Toolbar: FC = () => {
  return (
    <div className="absolute top-0 z-30 flex h-20 max-h-20 w-full flex-row items-center justify-start gap-5 bg-dark2 px-20">
      <FileUploadButton />
    </div>
  )
}
