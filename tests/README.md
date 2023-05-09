### How to test your code?
This system of tests uses ![pytest](https://docs.pytest.org/en/7.3.x/)
1. Always write tests in the `tests` folder
2. The name of your test script must start with `test_` (this is default and <b><i>can</i></b> be changed but for sake of brevity we will let it be)
3. Use the `import` statement to call your module
4. Within your script (z.B. test_createSVG.py), declare a function that starts with `test_*`
5. You may have as many functions you want but the ones that are to be tested should start
with `test_*`
6. Now call your module within the function by providing some dummy inputs
7. Remember these inputs can be intentionally wrong or right
8. Write `assert` statements to check if your input returns the correct result

### How to run your pytest?
>pytest tests/test_createSVG.py

### How to run all pytest?
>pytest tests/*