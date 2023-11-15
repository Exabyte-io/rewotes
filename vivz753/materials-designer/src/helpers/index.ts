import { Vector3 } from "three"

/**
 * Parses comma separated VECTORS (not values, as CSV usually stands for) into an array of 3d vectors
 * @param input text read from the user input
 * @returns an array containing vectors
 */
export const parseCSV = (input: string): Vector3[][] => {
  const lines = input.split("\n\n")
  const lines2 = lines.map((l) => l.split("\n"))
  const tokens = lines2.map((s) => s.map((s2) => s2.split(",")))
  const vectors = tokens.map((line) => line.map((l) => new Vector3(parseFloat(l[0]), parseFloat(l[1]), parseFloat(l[2]))))
  return vectors
}

/**
 * Parses .xyz files into 2 arrays, one for 3d vectors, one for elemental symbols
 * @param input text read from an .xyz file
 * @returns an array containing the vector & symbol arrays
 */
export const parseXYZ = (input: string): Array<string[] | Vector3[]> => {
  const lines = input.trim().split("\n")
  const tokens = lines.map((l) => {
    const line = l.split(" ")
    return line
  })
  const vectors = tokens.map(
    (line) => line && new Vector3(parseFloat(line[1]), parseFloat(line[2]), parseFloat(line[3]))
  )
  const symbols = tokens.map((line) => line && line[0])
  return [symbols, vectors]
}
