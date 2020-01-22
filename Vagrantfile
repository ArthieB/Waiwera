#-*- mode: ruby -*-
# vi: set ft=ruby :

VAGRANTFILE_API_VERSION = "2"
Vagrant.require_version ">= 2.2.4"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
 config.vm.box = "bento/debian-9"
 config.vm.synced_folder ".", "/vagrant/"
 config.vm.define :waiwera, primary: true do |waiwera|
   waiwera.vm.provision :shell, inline: "apt install python-pip -y && pip install -U ansible"
   waiwera.vm.provision "ansible_local" do |ansible|
     #ansible.limit = 'waiwera'
     ansible.provisioning_path = "/vagrant/install/ansible"
     ansible.verbose = "v"
     ansible.playbook = "vagrant.yml"
     # ansible.raw_arguments = [
     # ]
   end
   waiwera.vm.provision :shell, inline: "cd /vagrant ; python unit_tests.py", privileged: false
   waiwera.vm.hostname = "waiwera-debian"
 end
 config.vm.provider :virtualbox do |v|
#    v.gui = true
   v.memory = 4096
 end
 if Vagrant.has_plugin?("vagrant-cachier")
   config.cache.scope = :box
 end
end

