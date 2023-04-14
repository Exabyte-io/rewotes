const Meteor = {
    isClient: true,
    isServer: false,
    methods: () => {},
    call: jest.fn(),
    startup: jest.fn(),
    subscribe: jest.fn(),
};
  
export default Meteor