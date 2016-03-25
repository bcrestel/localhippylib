classdef HippyClient
    properties
        conn
    end
    
    methods
        function self = HippyClient
            self.conn = tcpip('localhost',1983,'Timeout', inf, 'OutputBufferSize',1024,'InputBufferSize',1024);
            fopen(self.conn);
            
        end
        
        function out = unspace(self,s)
            out = s(s ~= ' ');
        end
        
        function close(self)
            fwrite(self.conn,sprintf('%32s', 'Quit'));
            fclose(self.conn);
        end

        function computeMapPoint(self)
            fwrite(self.conn,sprintf('%32s', 'ComputeMapPoint'));
            val = self.unspace(char(fread(self.conn,32))');
            assert(~isempty(val));
            disp(val);
        end

        function kdim = KLE_GaussianPost(self)
             fwrite(self.conn,sprintf('%32s', 'KLE_GaussianPost'));
             kdim = fread(self.conn,1,'double');
             assert(~isempty(kdim));
             disp(kdim);
        end
        
        function val = negLogPost(self, eta)
            fwrite(self.conn,sprintf('%32s', 'NegLogPost'));
            val = self.send_array_receive_value(eta);
        end
        
        function val = negLogLikelihood(self, eta)
            fwrite(self.conn,sprintf('%32s', 'NegLogLikelihood'));
            val = self.send_array_receive_value(eta);
        end
        
        function val = negLogPrior(self, eta)
            fwrite(self.conn,sprintf('%32s', 'NegLogPrior'));
            val = self.send_array_receive_value(eta);
        end
        
        function val = negLogGaussianPost(self, eta)
            fwrite(self.conn,sprintf('%32s', 'NegLogGaussianPost'));
            val = self.send_array_receive_value(eta);
        end
        
        function echo(self,eta)
            fwrite(self.conn,sprintf('%32s', 'Echo'));
            self.send_array_receive_value(eta);
        end
        
        function val = send_array_receive_value(self, eta)
            fwrite(self.conn,sprintf('%30s',num2str(size(eta))));
            n_floats = numel(eta);
%            disp('Sending eta');
            n_sent = 0;
            while n_sent < n_floats
                n_tosend = min(128,n_floats-n_sent);
                fwrite(self.conn,eta((n_sent+1):n_sent+n_tosend),'double');
                n_sent = n_sent + n_tosend;
%                fprintf('%i/%i floats sent\n',n_sent,n_floats);    
            end
            val = fread(self.conn,1,'double');
        end
            

%         function put(self,name,val)
%             fwrite(self.conn,sprintf('%10s','put'));
%             fwrite(self.conn,sprintf('%10s',name));
%             fwrite(self.conn,sprintf('%30s',num2str(size(val))));
%             n_floats = numel(val);
%             n_sent = 0;
%             while n_sent < n_floats
%                 n_tosend = min(128,n_floats-n_sent);
%                 fwrite(self.conn,val((n_sent+1):n_sent+n_tosend),'double');
%                 n_sent = n_sent + n_tosend;
%                 fprintf('%i/%i floats sent\n',n_sent,n_floats);  
%             end
%         end
%         
%         function execute(self,stmt)
%             fwrite(self.conn,sprintf('%10s','exec'));
%             fwrite(self.conn,sprintf('%100s',stmt));
%         end
%         
%         function val = get(self,name)
%             fwrite(self.conn,sprintf('%10s','get'));
%             fwrite(self.conn,sprintf('%10s',name));
%             shape_str = char(fread(self.conn, 30,'char'))';
%             disp(shape_str);
%             shape = eval(shape_str);
%             val = zeros(shape);
%             n_floats = numel(val);
%             n_read = 0;
%             while n_read < n_floats
%                 n_toread = min(128,n_floats-n_read);
%                 val(n_read+1:n_read+n_toread) = fread(self.conn,n_toread,'double');
%                 n_read = n_read + n_toread;
%                 fprintf('%i/%i floats read\n',n_read,n_floats)
%             end
%             
%         end
        
    end
end
