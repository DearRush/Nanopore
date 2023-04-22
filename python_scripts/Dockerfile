# 基于的基础镜像

FROM python:3.9.6

# 将docker_test目录下的代码添加到镜像中的code文件夹（两个目录参数中间有空格分开）

COPY /public/python_scripts/ /code

# 设置code文件夹是工作目录

WORKDIR /code

# 安装支持

RUN pip install -r new.txt

#docker运行时即运行app.py文件

CMD ["python","/code/Identification.py",""]