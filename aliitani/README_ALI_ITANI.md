## NOTES - Implementation

I have implemented the basic arithmetic units as flowchart elements, where a user can run increment, decrement, multiply, or divide operations on a value.

I also implemented conditions where during a workflow run, conditional checks run against the value, that are the basic conditional operators,

- == equal to
- != not equal to
- '>' greater than
- '>=' greater than or equal to
- < less than
- <= less than or equal to

For the UI, I have implemented 2 viewer components the JSON and the flowchart side-by-side and made sure the JSON viewer that includes the JSON data structure updates reactively based on the flowchart content.

I also added test result execution workflow, to test if the workflow output value is as the asserted value as requested.

## TECH STACK

For the tech stack I have went with React, Vite, Tailwind CSS, and Reactflow as recommended.

- I used Reactflow to get the flowchart integrated with React.
- I used Tailwind CSS to get a more styled UI.
- I ended up not using meteor as suggested, since this app is solely a ui application without any backend logic required. Based on that, I ended up going with React+Vite for the frontend app.

## TO RUN THE Application

Open terminal and make sure you are in the project app directory, and run the following command after running npm or yarn to install npm dependencies.

> $ npm run dev
