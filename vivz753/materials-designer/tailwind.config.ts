import type { Config } from "tailwindcss"

const config: Config = {
  content: [
    "./src/pages/**/*.{js,ts,jsx,tsx,mdx}",
    "./src/components/**/*.{js,ts,jsx,tsx,mdx}",
    "./src/app/**/*.{js,ts,jsx,tsx,mdx}",
  ],
  theme: {
    extend: {
      fontFamily: {
        pixelboy: ["Pixelboy"],
        oldschool: ["OldSchoolAdventures"],
        mozart: ["Mozart"],
        neoeuler: ["NeoEuler"],
      },
      colors: {
        dark3: "#121616",
        dark2: "#283333",
        dark1: "#3e5050",
        accent: "#90c641",
        light: "#f4f5f6",
      },
      backgroundImage: {
        "gradient-radial": "radial-gradient(var(--tw-gradient-stops))",
        "gradient-conic": "conic-gradient(from 180deg at 50% 50%, var(--tw-gradient-stops))",
      },
    },
  },
  plugins: [],
}
export default config
