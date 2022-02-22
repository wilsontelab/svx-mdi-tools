---
published: false
---

## Create a Tier 3 AMI for public MDI tool servers

Follow these instructions to create a Tier 3 AMI to deliver
the tools, especially the Stage 2 apps, in this suite.

See general instructions for Tier 3 AMI creation here:

- <https://github.com/MiDataInt/mdi-aws-ami>

Please note: do NOT edit server.sh yet, that is done per 
server instance, not in the AMI.

### Initialize parent instance

Launch a new AWS EC2 instance from an appropriate Tier 2 empty AMI,
with the following properties:

- **instance type** = t3.xlarge
- **storage** = 20 GB SSD

### Install this tool suite

Add GIT_USER/SUITE_NAME to the suites list:

```bash
server edit suites.yml
```

Then:

```bash
server build
server install
```

### Download/install resources (optional)

```bash
cd /srv/mdi/resource-scripts
git clone https://github.com/GIT_USER/RESOURCES_NAME.git
server resource RESOURCES_NAME/SCRIPT_NAME.sh
```

### Save the image

Continue following the general instructions linked above 
to clean and create the public AWS AMI.
