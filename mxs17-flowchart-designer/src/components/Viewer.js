import { EDITOR_H, PALETTE_W } from "../constants";

export default function Viewer({ data }) {
  return (
    <div style={{ maxHeight: EDITOR_H, width: PALETTE_W * 3 }}>
      <pre>{data}</pre>
    </div>
  );
}
