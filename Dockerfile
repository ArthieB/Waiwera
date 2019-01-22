FROM debian:9

ARG waiwera_user
ARG waiwera_pwd

ENV PETSC_DIR=/opt/app/waiwera/external/PETSc
ENV PETSC_ARCH=arch-linux2-c-debug
ENV LD_LIBRARY_PATH="/opt/app/lib"
ENV PYTHONPATH=/"/opt/app/PyTOUGH:/opt/app/credo2:/opt/app/supermodels-test/utils"
ENV PATH="$PATH:/opt/app/bin:/opt/app/supermodels-test/dist"

RUN apt update && \
    apt install -y python3-minimal python3-pip && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*  /root/.cache/*

RUN pip3 install pip
RUN pip3 install ansible

ADD ansible /srv/ansible
WORKDIR /srv/

RUN /usr/local/bin/ansible-playbook -c local ansible/site.yml -v -e "waiwera_user=${waiwera_user}" -e "waiwera_pwd=${waiwera_pwd}"

WORKDIR /opt/app/supermodels-test

RUN pip3 uninstall ansible -y && \
    apt-get autoremove -y && apt-get autoclean -y && apt-get clean -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cache/* /srv/ansible /usr/share/doc/*

RUN for dep in $(pip show ansible | grep Requires | sed 's/Requires: //g; s/,//g'); do pip uninstall -y $dep; done
