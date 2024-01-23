import { useEffect, useRef } from "react";

export default function BlockInput({ top, left, width, height, submitValue }) {
  const inputRef = useRef();

  useEffect(() => {
    inputRef.current.focus();
  }, []);

  return (
    <input
      ref={inputRef}
      style={{ top, left, width, height }}
      className="num-input"
      onBlur={(e) => submitValue(e.target.value)}
    />
  );
}
