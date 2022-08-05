print('Imported my module...')

Test = 'Test string'


def find_index(to_search, target):
    '''Find the index for a value in a sequence'''
    for i, value in enumerate(to_search):
        if value == target:
            return i
    return -1
