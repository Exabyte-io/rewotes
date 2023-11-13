import styled from 'styled-components';

export const DesignerContainer = styled.div`
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
`;

export const JSONContainer = styled.pre`
  margin-left: 20px;
  margin-top: 0;
  width: 20vw;
  max-height: 60vh;
  overflow: auto;
  border: 1px solid #ccc;
  padding: 10px;
`;

export const FlowchartElement = styled.div<{ $isDisabled: boolean }>`
  max-width: 120px;
  height: 60px;
  border: 1px solid black;
  margin: 10px 4px;
  display: flex;
  align-items: center;
  justify-content: center;
  background-color: #00435a;
  opacity: ${({ $isDisabled }) => ($isDisabled ? 1 : 0.2)};
  color: white;
  border-radius: 12px;
  padding: 4px;
`;

export const StyledAside = styled.aside`
  display: flex;
  width: 100%;
  overflow-x: auto;
`;

export const SaveButton = styled.button`
  margin-left: 4%;
  border-radius: 12px;
  borber: 1px solid #cc5500;
  color: black;
  font-weight: 600;
  background-color: ${(props) => (props.disabled ? 'rgba(204, 85, 0, 0.2)' : '#CC5500')};
`;

export const ResetButton = styled(SaveButton)`
  margin: 10px 0 10px auto;
  width: 120px;
  height: 60px;
`;

export const CanvasWrapper = styled.div`
  height: 62vh;
  width: 60vw;
  background-color: rgba(0, 67, 90, 0.9);
  border: 1px solid #ccc;
`;
