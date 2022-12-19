import { useState } from 'react';
import { useSnapshot } from 'valtio';
import { JsonDataState } from '../state/AppState';

/***
 * This component holds the test case evaluation, that you can run against the workflow executed result.
 * Once the workflow is executed and output value is generated based on the result, then the assertion is evaluated.
 */
const TestExecution = () => {
	const [expectedValue, setExpectedValue] = useState(0);
	const json = useSnapshot(JsonDataState);
	const result =
		expectedValue === json.outputValue ? (
			<span className='text-green-500'>PASSED</span>
		) : (
			<span className='text-red-500'>FAILED</span>
		);

	return (
		<div className='flex flex-col flex-nowrap gap-2 w-full bg-slate-600 p-2 border-l-blue-600 border-l-4'>
			<p>
				<span className='font-bold italic'>Test Execution Result</span>
			</p>
			<hr />
			<p>Test the workflow outputValue against an expectedValue</p>
			<p>The outputValue is, {json.outputValue}</p>
			<label htmlFor='test-case'>
				Expected Value:
				<input
					className='text-slate-700 mx-2 text-right'
					type='number'
					id='test-case'
					value={expectedValue}
					onChange={(e) => setExpectedValue(Number(e.target.value))}
				/>
			</label>
			<p>
				Assertion result of the expected value and actual outputValue:{' '}
				{json.isExecutionDone && result}
			</p>
		</div>
	);
};

export default TestExecution;
