FROM xnat-env-bootstrap:brats

#just overwrite changes to pymipl sources if any.

# Bootstrap entrypoint
COPY ./xnat_env_bootstrap/main.sh /opt/entrypoint/
COPY *.py /opt/packages/pymipl/
COPY xnat_workflow /opt/packages/pymipl/xnat_workflow

RUN chmod 0755 /opt/entrypoint/main.sh && chmod 777 /workdir && chmod 777 /opt/packages/user

#CMD ["/opt/bootstrap/bootstrap.sh"]
