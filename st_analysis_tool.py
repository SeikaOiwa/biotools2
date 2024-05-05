import streamlit as st
from PIL import Image
import os

base_f_path_ = os.getcwd()
base_f_path = f'{base_f_path_}/File'

def main():
    st.image(f'{base_f_path}/main.png')   
    st.image(f'{base_f_path}/main2.png',use_column_width=True)  
    
if __name__ == "__main__":
    main()