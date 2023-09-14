import { ChangeEvent, FC, MouseEventHandler, useRef } from "react"

export const FileUploadButton: FC<{ handleFile?: (file: any) => void }> = ({ handleFile }) => {
  const fileInput = useRef<HTMLInputElement>(null)

  const handleClick: MouseEventHandler<HTMLButtonElement> = (event) => {
    fileInput.current?.click()
    console.log()
  }

  const handleUpload = (event: ChangeEvent<HTMLInputElement>) => {
    const fileUploaded = event.target.files && event.target.files[0]
    handleFile && handleFile(fileUploaded)
    // TODO: handle the file somehow
    console.log(fileUploaded)
  }

  return (
    <button
      className="flex items-center justify-center rounded-md bg-accent px-2 py-0.5 font-mozart text-xl uppercase tracking-widest text-light focus:outline-light"
      onClick={handleClick}
    >
      Import
      {/* TODO: accept specific file type only */}
      <input ref={fileInput} onChange={handleUpload} type="file" accept=".xyz,.poscar,.txt" className="hidden" />
    </button>
  )
}
